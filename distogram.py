"""
This script computes a distogramm for an MD trajectory 
with weights.
"""
import numpy as np
import scipy.spatial
import mdtraj as md
import numba as nb
import matplotlib.pyplot as plt
import tqdm # not necessary, it's jsut for progress visualization of a for loop




@nb.jit
def increment_counts(counts, dist_bin_int, weight):
    n_ca, _ = dist_bin_int.shape
    for i in range(n_ca):
        for j in range(n_ca):
            b = dist_bin_int[i, j]
            counts[i, j, b] += weight
    return counts



########## input params and input data paths ###################

topPath = "topology.pdb"
trajPath = "trajectory.dcd"
weightPath =  'conformer_weights.dat'


n_dist = 32 # number of distance bins
range_dist = 0, 150 # distance range

#################################################################

# load input data

weights = np.loadtxt(weightPath)[:, 1]

t = md.load(trajPath, top=topPath)
n_frames = len(t)
top = t.top
ca_idx = top.select("protein and name CA")
n_ca = len(ca_idx)




axis_dist = np.linspace(*range_dist, n_dist)
axis_spacing = np.diff(axis_dist)[0] # bin width

counts = np.zeros([n_ca, n_ca, n_dist], dtype=np.float32)
counts += 1e-12 # add a small offset to avoid NaNs in later scoring

# distogram calculation
for i_frame in tqdm.tqdm(range(n_frames)):
    weight, frame = weights[i_frame], t[i_frame]
    xyz_ca = frame.xyz[0, ca_idx,: ]
    dist = scipy.spatial.distance_matrix(xyz_ca, xyz_ca) * 10.0 # units in mdtraj in nm
    # find bin 
    dist_bin = dist / axis_spacing
    dist_bin_int = dist_bin.astype(dtype=np.int32)
    counts = increment_counts(counts, dist_bin_int, weight)



fig, ax = plt.subplots() 
sel_bin = 1 # distogram is multidim. array hard to visualize -> instead you can print 2D map at selected distance bin
im = plt.imshow(counts[:, :, sel_bin])
ax.set_xlabel(r'$C_{\alpha}\ $residue number')
ax.set_ylabel(r'$C_{\alpha}\ $residue number')

ax.set_title("distance bin = %.2f Å" % (axis_dist[sel_bin]))
plt.colorbar(im,orientation='vertical').set_label(label="weighted distance occupancy",size=10)


plt.savefig('occupancy_distance_bin_%.2fA.png'% (axis_dist[sel_bin]),dpi=300)    
plt.close(fig)

# save distogramm
np.save("distogram.npy", counts)
np.save("distogram_axis.npy", axis_dist)


# Compute mean and SD and skewness of distogram
means = np.zeros((n_ca, n_ca), dtype=np.float32)
sds = np.zeros((n_ca, n_ca), dtype=np.float32)
skws = np.zeros((n_ca, n_ca), dtype=np.float32)
axis_dist2 = axis_dist**2.0
axis_dist3 = axis_dist**3.0
for i in range(n_ca):
    for j in range(n_ca):
        s = np.clip(np.sum(counts[i, j, :]), 1e-12, 1.0) # "normalize" histogram to 1
        mean = np.dot(counts[i, j, :], axis_dist) / s 
        mean2 = np.dot(counts[i, j, :], axis_dist2) / s
        mean3 = np.dot(counts[i, j, :], axis_dist3) / s
        d2 = mean2 - mean**2.0
        sd = np.sqrt(abs(d2))
        skw = (mean3 -  3 * mean * sd**2.0 - mean**3.0) # / (sd**3.0) # removed normalization (numerical unstable)
        means[i, j] = mean
        sds[i, j] = sd
        skws[i, j] = skw

r = np.tanh(abs(means / np.max(means)))
g = np.tanh(abs(sds / np.max(sds)))
b = np.tanh(abs(skws / np.max(skws)))

rgb = np.dstack((r,g,b))

fig, ax = plt.subplots() 
im2 = plt.imshow(rgb) 
ax.set_xlabel(r'$C_{\alpha}\ $residue number')
ax.set_ylabel(r'$C_{\alpha}\ $residue number')
plt.savefig('mu_sds_skew_of_distogram.png',dpi=300)    
plt.close(fig)
