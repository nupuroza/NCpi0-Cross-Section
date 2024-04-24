// Classes to convert between THnD and TMatrix for 1 and 2 dimensions

// Create TMatrixD from TH2D
// `invert_matrix` allows for the vertical/horizontal axis convention of the input TH2D to be reversed
// `include_{under,over}flow_{X,Y}` enables the user to include or exclude {under,over}flow bins from conversion
TMatrixD TH2DtoTMatrixD(const TH2D& input_hist, bool invert_matrix, bool include_underflow_X, bool include_overflow_X, bool include_underflow_Y, bool include_overflow_Y)
{
  // Number of bins are increased by one unit for each of the underflow and overflow that is included
  Int_t nRows = input_hist.GetNbinsX();
  Int_t nColumns = input_hist.GetNbinsY();
  
  if(include_underflow_X) nRows = nRows+1;
  if(include_overflow_X)  nRows = nRows+1;

  if(include_underflow_Y) nColumns = nColumns+1;
  if(include_overflow_Y)  nColumns = nColumns+1;
   
  TMatrixD output_matrix(nRows,nColumns);
  // The following loops start at zero only if we're including underflow
  Int_t offset_X = include_underflow_X ? 0 : 1;  
  Int_t offset_Y = include_underflow_Y ? 0 : 1;
  for(Int_t i=0; i<nRows; i++)
  {
    for(Int_t j=0; j<nColumns; j++)
    {
      if(invert_matrix) output_matrix(i, j) = input_hist.GetBinContent(j+offset_Y, i+offset_X);
      else              output_matrix(i, j) = input_hist.GetBinContent(i+offset_X, j+offset_Y);
    }
  }
  return output_matrix;
}

// Create TMatrixD from TH1D
// `include_{under,over}flow` enables the user to include or exclude {under,over}flow bins from conversion
TMatrixD TH1DtoTMatrixD(const TH1D& input_hist, bool include_underflow, bool include_overflow)
{
  // Number of bins are increased by one unit for each of the underflow and overflow that is included
  Int_t nRows = input_hist.GetNbinsX();

  if(include_underflow) nRows = nRows+1;
  if(include_overflow)  nRows = nRows+1;

  TMatrixD output_matrix(nRows,1);
  // The following loop starts at zero only if we're including underflow
  Int_t offset = include_underflow ? 0 : 1;  
  for(Int_t i=0; i<nRows; i++)
    {
      output_matrix(i, 0) = input_hist.GetBinContent(i+offset);
    }
  return output_matrix;
}

// Create TH2D from TMatrixD
// `reference_hist` should have the desired binning, but will not itself be filled
// DISCLAIMER: As implemented this evidently assumes you want to convert a square matrix
// It is conceivably the case that you might want to make reference_hist a TH2D in the future
TH2D TMatrixDtoTH2D(const TMatrixD& input_mat, const TH1D& reference_hist)
{
  Int_t nBins = reference_hist.GetNbinsX();
  Double_t binEdges[nBins+1];
  for(int i=0; i<nBins+1; i++){
      binEdges[i] = reference_hist.GetBinLowEdge(i+1);
  }
  TH2D output_hist("","",nBins,binEdges,nBins,binEdges);
  for(Int_t i=0; i<input_mat.GetNrows(); i++)
    {
      for(Int_t j=0; j<input_mat.GetNcols(); j++)
      {
	      output_hist.SetBinContent(i+1, j+1, input_mat(j, i));
      }
    }
  return output_hist;
}

// Create TH1D from TMatrixD
// `reference_hist` should have the desired binning, but will not itself be filled
TH1D TMatrixDtoTH1D(const TMatrixD& input_mat, const TH1D& reference_hist)
{
  Int_t nBins = reference_hist.GetNbinsX();
  Double_t binEdges[nBins+1];
  for(int i=0; i<nBins+1; i++){
      binEdges[i] = reference_hist.GetBinLowEdge(i+1);
  }
  TH1D output_hist("","",nBins,binEdges);
  for(Int_t i=0; i<input_mat.GetNrows(); i++)
    {
      output_hist.SetBinContent(i+1, input_mat.operator()(i, 0));
    }
  return output_hist;
}

