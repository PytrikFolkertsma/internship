import loompy
import glob
import velocyto as vcy
import numpy as np
from sklearn.manifold import TSNE


print('loading data')
vlm = vcy.load_velocyto_hdf5('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/velocyto/180831.hdf5')

print(len(vlm.ca['CellID']), 'cells')
print(len(vlm.ra['Gene']), 'genes')

#print('filtering cells')
#vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))

print('filtering genes')
#vlm.score_cv_vs_mean(3000, plot=False, max_expr_avg=35)
vlm.score_cv_vs_mean(3000, plot=False)
vlm.filter_genes(by_cv_vs_mean=True)

print(len(vlm.ca['CellID']), 'cells')
print(len(vlm.ra['Gene']), 'genes')

print('normalizing data matrices')
vlm._normalize_S(relative_size=vlm.S.sum(0), target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0), target_size=vlm.U.sum(0).mean())

print('running pca')
vlm.perform_PCA()

print('knn smoothing')
#vlm.knn_imputation(n_pca_dims=15, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=20)
vlm.knn_imputation(n_pca_dims=20, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=20)

#print('fit gammas')
vlm.fit_gammas()

print('calculate velocity')
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)

print('running tsne')
bh_tsne = TSNE(random_state=1)
vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :20])

print('projection of velocity onto embeddings')
vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1, n_neighbors=3500, knn_random=True, sampled_fraction=0.5)

print('calculate embedding shift')
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

print('calculate grid arrows')
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=100)

print('saving hdf5')
vlm.to_hdf5('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/velocyto/180831_tsne1_velocity.hdf5')






