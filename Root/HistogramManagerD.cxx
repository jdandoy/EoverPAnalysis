/******************************************
 *
 * Base class used to book a set of histograms.
 * Many sets can be created and compared for
 * efficiecny studies or population studies.
 *
 * G. Facini
 * Sep 9 08:58:12 CEST 2014
 * G. Stark
 * Jan 20 10:29 CEST 2015
 *
 ******************************************/


#include "EoverPAnalysis/HistogramManagerD.h"

/* constructors and destructors */
HistogramManagerD::HistogramManagerD(std::string name, std::string detailStr):
  m_name(name),
  m_detailStr(detailStr)
{

  // if last character of name is a alphanumeric add a / so that
  // in the output file, a TDirectory is created with the histograms inside
  if( isalnum( m_name.back() ) && !ispunct( m_name.back() ) ) {
    m_name += "/";
    //Info("HistogramManagerD()", "Adding slash to put hists in TDirectories: %s",m_name.c_str());
  }

}

HistogramManagerD::~HistogramManagerD() {}

/* Main book() functions for 1D, 2D, 3D histograms */
TH1D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, double xlow, double xhigh)
{
  TH1D* tmp = new TH1D( (name + title).c_str(), title.c_str(), xbins, xlow, xhigh);
  SetLabel(tmp, xlabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, double xlow, double xhigh,
                             std::string ylabel, int ybins, double ylow, double yhigh)
{
  TH2D* tmp = new TH2D( (name + title).c_str(), title.c_str(), xbins, xlow, xhigh, ybins, ylow, yhigh);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH3D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, double xlow, double xhigh,
                             std::string ylabel, int ybins, double ylow, double yhigh,
                             std::string zlabel, int zbins, double zlow, double zhigh)
{
  TH3D* tmp = new TH3D( (name + title).c_str(), title.c_str(), xbins, xlow, xhigh, ybins, ylow, yhigh, zbins, zlow, zhigh);
  SetLabel(tmp, xlabel, ylabel, zlabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TProfile* HistogramManagerD::book(std::string name, std::string title,
				 std::string xlabel, int xbins, double xlow, double xhigh,
				 std::string ylabel, double ylow, double yhigh, 
				 std::string option)
{
  TProfile* tmp = new TProfile( (name + title).c_str(), title.c_str(), xbins, xlow, xhigh, ylow, yhigh, option.c_str());
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TProfile* HistogramManagerD::book(std::string name, std::string title,
				 std::string xlabel, int xbins, const Double_t *xbinsArr,
				 std::string ylabel, double ylow, double yhigh, 
				 std::string option)
{
  TProfile* tmp = new TProfile( (name + title).c_str(), title.c_str(), xbins, xbinsArr, ylow, yhigh, option.c_str());
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

/////// Variable Binned Histograms ///////
TH1D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, const Double_t* xbinArr)
{
  TH1D* tmp = new TH1D( (name + title).c_str(), title.c_str(), xbins, xbinArr);
  SetLabel(tmp, xlabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, const Double_t* xbinArr,
                             std::string ylabel, int ybins, double ylow, double yhigh)
{
  TH2D* tmp = new TH2D( (name + title).c_str(), title.c_str(), xbins, xbinArr, ybins, ylow, yhigh);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, double xlow, double xhigh,
                             std::string ylabel, int ybins, const Double_t* ybinArr)
{
  TH2D* tmp = new TH2D( (name + title).c_str(), title.c_str(), xbins, xlow, xhigh, ybins, ybinArr);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, const Double_t* xbinArr,
                             std::string ylabel, int ybins, const Double_t* ybinArr)
{
  TH2D* tmp = new TH2D( (name + title).c_str(), title.c_str(), xbins, xbinArr, ybins, ybinArr);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH3D* HistogramManagerD::book(std::string name, std::string title,
                             std::string xlabel, int xbins, const Double_t* xbinArr,
                             std::string ylabel, int ybins, const Double_t* ybinArr,
                             std::string zlabel, int zbins, const Double_t* zbinArr)
{
  TH3D* tmp = new TH3D( (name + title).c_str(), title.c_str(), xbins, xbinArr, ybins, ybinArr, zbins, zbinArr);
  SetLabel(tmp, xlabel, ylabel, zlabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

/* Helper functions */
void HistogramManagerD::Sumw2(TH1* hist, bool flag /*=true*/) {
  hist->Sumw2(flag);
}

void HistogramManagerD::record(TH1* hist) {
  m_allHists.push_back( hist );
}

void HistogramManagerD::record(EL::Worker* wk) {
  for( auto hist : m_allHists ){
    wk->addOutput(hist);
  }
}

void HistogramManagerD::SetLabel(TH1* hist, std::string xlabel)
{
  hist->GetXaxis()->SetTitle(xlabel.c_str());
}

void HistogramManagerD::SetLabel(TH1* hist, std::string xlabel, std::string ylabel)
{
  hist->GetYaxis()->SetTitle(ylabel.c_str());
  this->SetLabel(hist, xlabel);
}

void HistogramManagerD::SetLabel(TH1* hist, std::string xlabel, std::string ylabel, std::string zlabel)
{
  hist->GetZaxis()->SetTitle(zlabel.c_str());
  this->SetLabel(hist, xlabel, ylabel);
}
