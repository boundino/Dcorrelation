#ifndef _THGRSTYLE_H_
#define _THGRSTYLE_H_

/* thgrstyle.h */

namespace xjjroot
{
  class thgrstyle
  {
  public:
    thgrstyle(Color_t mcolor=-1, Style_t mstyle=-1, Size_t msize=-1, Color_t lcolor=-1, Style_t lstyle=-1, Width_t lwidth=-1, Color_t fcolor=-1, Float_t falpha=-1, Style_t fstyle=-1, TString opt="")
      : m_mcolor(mcolor)
      , m_mstyle(mstyle)
      , m_msize(msize)
      , m_lcolor(lcolor)
      , m_lstyle(lstyle)
      , m_lwidth(lwidth)
      , m_fcolor(fcolor)
      , m_falpha(falpha)
      , m_fstyle(fstyle)
      , m_opt(opt) 
    {}        
    ~thgrstyle() {};
    Color_t GetMarkerColor() const {return m_mcolor;}
    Style_t GetMarkerStyle() const {return m_mstyle;}
    Size_t GetMarkerSize() const {return m_msize;}
    Color_t GetLineColor() const {return m_lcolor;}
    Style_t GetLineStyle() const {return m_lstyle;}
    Width_t GetLineWidth() const {return m_lwidth;}
    Color_t GetFillColor() const {return m_fcolor;}
    Float_t GetFillAlpha() const {return m_falpha;}
    Style_t GetFillStyle() const {return m_fstyle;}
    TString GetOption() const {return m_opt;}
  private:
    Color_t m_mcolor;
    Style_t m_mstyle;
    Size_t m_msize;
    Color_t m_lcolor;
    Style_t m_lstyle;
    Width_t m_lwidth;
    Color_t m_fcolor;
    Float_t m_falpha;
    Style_t m_fstyle;
    TString m_opt;
  };
  
  template <class T> void setthgrstyle(T* h, const thgrstyle &sty);

}

template <class T> 
void xjjroot::setthgrstyle(T* h, const thgrstyle &sty)
{
  if(sty.GetMarkerColor()>=0) h->SetMarkerColor(sty.GetMarkerColor());
  if(sty.GetMarkerStyle()>=0) h->SetMarkerStyle(sty.GetMarkerStyle());
  if(sty.GetMarkerSize()>=0)  h->SetMarkerSize(sty.GetMarkerSize());
  if(sty.GetLineColor()>=0) h->SetLineColor(sty.GetLineColor());
  if(sty.GetLineStyle()>=0) h->SetLineStyle(sty.GetLineStyle());
  if(sty.GetLineWidth()>=0) h->SetLineWidth(sty.GetLineWidth());
  if(sty.GetFillColor()>=0) h->SetFillColor(sty.GetFillColor());
  if(sty.GetFillAlpha()>=0) h->SetFillColorAlpha(sty.GetFillColor(), sty.GetFillAlpha());
  if(sty.GetFillStyle()>=0) h->SetFillStyle(sty.GetFillStyle());
}

#endif

