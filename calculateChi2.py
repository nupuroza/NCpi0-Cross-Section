### Script to calculate chi squared values from.
### Uses the Combined Neyman-Pearson method taken from https://arxiv.org/abs/1903.07185

import ROOT
import sys

# mc_hist is a TH1 with the generator prediction
# data_hist is a TH1 with the data distribution being compared against mc_hist
# cov_matrix is a TMatrix with the data covariance matrix
# cov_matrix_includes_stats is a boolean indicating whether cov_matrix already includes statistical uncertainties. If not, it will be calculated. Can also be set to True if you don't want to include statistical uncertainties.
# test comment
def calculateChi2(mc_hist, data_hist, cov_matrix, cov_matrix_includes_stats):
    diff_hist = data_hist.Clone("diff_hist")
    diff_hist.Add(mc_hist, -1)
    # Overflow bins currently excluded from analysis. Will likely be changed at a later date.
    nbins = diff_hist.GetNbinsX()
    if nbins + 1 != cov_matrix.GetNrows() or nbins + 1 != cov_matrix.GetNcols():
        print("nbins + 1: " + str(nbins + 1))
        print("cov_matrix_rows: " + str(cov_matrix.GetNrows()))
        print("cov_matrix_cols: " + str(cov_matrix.GetNcols()))
        print("ERROR: Histogram binning does not match covariance matrix in calculateChi2(). Exiting!")
        sys.exit(1)
    # If statistical uncertainties aren't already included, calculate them using Eq. 19 of paper.
    if not cov_matrix_includes_stats:
        for i in range(nbins + 1):
            cov_matrix[i][i] += 3/(1/data_hist.GetBinContent(i+1) + 2/mc_hist.GetBinContent(i+1))
    inv_cov_matrix = ROOT.TMatrixD(ROOT.TMatrixD.kInverted, cov_matrix)
    diff_mat = ROOT.TMatrixD(nbins + 1, 1)
    for i in range(nbins + 1):
        diff_mat[i][0] = diff_hist.GetBinContent(i + 1)
    # Calculate chi2 using Eq. 18 of paper.
    chi2_temp = ROOT.TMatrixD(diff_mat, ROOT.TMatrixD.kTransposeMult, inv_cov_matrix)
    chi2 = ROOT.TMatrixD(chi2_temp, ROOT.TMatrixD.kMult, diff_mat)
    return chi2(0, 0), cov_matrix
