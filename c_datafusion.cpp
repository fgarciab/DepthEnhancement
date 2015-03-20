/*
* Author: Frederic GARCIA BECERRO
* Email: frederic.garcia.becerro@gmail.com
* Website: http://www.frederic-garcia-becerro.com
*/

#include "c_datafusion.h"

//#define MLF_ANALYSE_TIME_CONSUMPTION 1

c_DataFusion::c_DataFusion(s16 p_sImageWidth, s16 p_sImageHeight, s16 p_sNr_shift)
{            
    m_DataFusion.m_sInitOk = 0;
    m_DataFusion.m_w_original = p_sImageWidth;
    m_DataFusion.m_h_original = p_sImageHeight;
    m_DataFusion.m_w = (p_sImageWidth>>p_sNr_shift);
    m_DataFusion.m_h = (p_sImageHeight>>p_sNr_shift);
    m_DataFusion.m_nr_shift = p_sNr_shift;
    u08 l_ucShiftUnit = (1<<m_DataFusion.m_nr_shift);  // 1..16
    // Fill Dtt table
    for (s16 yy= 0; yy< l_ucShiftUnit; yy++)
    for (s16 xx= 0; xx< l_ucShiftUnit; xx++)
    {
        m_DataFusion.m_saDxx[yy][xx][0] = (l_ucShiftUnit-xx) * (l_ucShiftUnit-yy);
        m_DataFusion.m_saDxx[yy][xx][1] = xx * (l_ucShiftUnit-yy);
        m_DataFusion.m_saDxx[yy][xx][2] = yy * (l_ucShiftUnit-xx);
        m_DataFusion.m_saDxx[yy][xx][3] = yy * xx;
    }        
    m_DataFusion.m_Set.m_fSigma_spatial = 10.0;
    m_DataFusion.m_Set.m_fSigma_range = 10.0;
    m_DataFusion.m_Set.m_fSigma_Q = 10.0;
    m_DataFusion.m_Set.m_sScaleFactor = 1<<m_DataFusion.m_nr_shift;

    // Gaussian LUT for the fusion filters
    FillColorWeightedTable(m_DataFusion.m_ucaTable, m_DataFusion.m_Set.m_fSigma_range, MLF_APP_2D_RANGE);
    // Gaussian LUT for the fusion filters (BF)
    FillColorWeightedTable(m_DataFusion.m_ucaTableBF, m_DataFusion.m_Set.m_fSigma_range, MLF_APP_INTENSITY_RANGE);
    // Initialize the Gaussian filter
    Gaussian_Recursive_Order0_Init(m_DataFusion.m_Set.m_fSigma_spatial/(f32)(m_DataFusion.m_Set.m_sScaleFactor), m_DataFusion.m_h, m_DataFusion.m_w, &m_DataFusion.m_StGaussian);
    // Initialize memory handling
    m_StMem.m_lMemUsed = 0;
    m_StMem.m_pMemPointer = m_cMem;
    m_StMem.m_lMemSize = sizeof(m_cMem);
    s32 l_lRet;

    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pdaImgBox,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_J_BF_DS,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(u16));
//FGa, RGB-D Filter
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_RGB,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08)*3);
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_GRAY,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_R,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_G,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_B,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_GRAY_DS,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_R_DS,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_G_DS,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_I_B_DS,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_D,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u16));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_D_DS,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(u16));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Q,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Q_DS,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_BetaValue,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_BetaChannel,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u08));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_J,m_DataFusion.m_w_original*m_DataFusion.m_h_original*sizeof(u16));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_R,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_R,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_G,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_G,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_B,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_B,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_GRAY,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_GRAY,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Wk_R,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Wk_G,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Wk_B,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Wk_GRAY,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
#ifdef	MLF_BackGround_Counter
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_R_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_R_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_G_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_G_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_B_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_B_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level0_GRAY_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
    l_lRet = Alloc_Mem(&m_StMem,(char**)&m_DataFusion.m_pt_Jk_level1_GRAY_BG,m_DataFusion.m_w*m_DataFusion.m_h*sizeof(s32));
#endif
    m_DataFusion.m_sInitOk = 1;
}

c_DataFusion::~c_DataFusion()
{
    // Nothing to be done yet...
}

s32 c_DataFusion::Alloc_Mem(T_Memory_St* p_pSt,char **p_ppVoid, s32 p_lMemSize)
{
    if (p_pSt->m_lMemUsed+p_lMemSize > p_pSt->m_lMemSize)
    {
        //Not enough memory available
        std::cerr << "c_DataFusion::Alloc_Mem, Not enought memory." << std::endl << "MemUsed: " << p_pSt->m_lMemUsed << " MemElt: " << p_lMemSize << " MemSize: " << p_pSt->m_lMemSize << std::endl;
        return -1;
    }
    *p_ppVoid = p_pSt->m_pMemPointer+p_pSt->m_lMemUsed;
    p_pSt->m_lMemUsed += p_lMemSize;

    return p_pSt->m_lMemUsed;
}

void c_DataFusion::FillColorWeightedTable(u08* p_pucTable, f32 p_pfSigma_range, s16 len)
{
    s16 y;    
    f32 l_fSigma_rangeSq = (f32)2*p_pfSigma_range*p_pfSigma_range;
    for (y=0;y<len;y++)
    {
        (*p_pucTable++) = (u08)(exp((f32)(-(y*y))/(l_fSigma_rangeSq))*255);
        /* Gaussian value between 0..255*/
    }
}

void c_DataFusion::Gaussian_Recursive_Order0_Init(f32 sigma, s16 h, s16 w, T_Gaussian_Rec_St* p_StOut)
{
    f32		alpha ;
    f32		ema ;
    f32		ema2;
    f32		nsigma;

    if (sigma < 0.1f)
        nsigma = 0.1f;
    else
        nsigma = sigma;
    alpha = 1.695f / nsigma;
    ema = exp(-alpha);
    ema2 = exp(-2*alpha);
    p_StOut->b1 = -2*ema;
    p_StOut->b2 = ema2;
    {
        const f32 k = (1-ema)*(1-ema)/(1+2*alpha*ema-ema2);
        p_StOut->a0 = k;
        p_StOut->a1 = k*(alpha-1)*ema;
        p_StOut->a2 = k*(alpha+1)*ema;
        p_StOut->a3 = -k*ema2;
    }

    p_StOut->coefp = (p_StOut->a0+p_StOut->a1)/(1+p_StOut->b1+p_StOut->b2);
    p_StOut->coefn = (p_StOut->a2+p_StOut->a3)/(1+p_StOut->b1+p_StOut->b2);
}

void c_DataFusion::Gaussian_Recursive_Order0(const T_Gaussian_Rec_St* p_StSet, s32 *p_pdaImg, s32 *p_pdaImgTemp, s16 h, s16 w)
{
    Gaussian_Recursive_x(p_pdaImgTemp,p_pdaImg,w,h,
                p_StSet->a0,p_StSet->a1,p_StSet->a2,p_StSet->a3,p_StSet->b1,p_StSet->b2,p_StSet->coefp,p_StSet->coefn);
    Gaussian_Recursive_y(p_pdaImg,p_pdaImgTemp,w,h,
                p_StSet->a0,p_StSet->a1,p_StSet->a2,p_StSet->a3,p_StSet->b1,p_StSet->b2,p_StSet->coefp,p_StSet->coefn);
}

void c_DataFusion::Gaussian_Recursive_x(s32 *p_pdImgOut,const s32 *p_pdImgIn, s16 w, s16 h,
            const f32 a0, const f32 p_dA1, const f32 a2, const f32 a3,
            const f32 b1, const f32 b2, const f32 coefp, const f32 coefn)
{
    const s32 *l_pdImgIn;
    s32 *l_pus_J;
    {
        s16 x;
        s16 y;
        s32 xp;  // previous input
        f32 yp;  // previous output
        f32 yb;  // previous output by 2
        l_pdImgIn =  p_pdImgIn;
        l_pus_J =  p_pdImgOut;

        for (y=0;y<h;y++)
        {
            x= w;
            xp = *l_pdImgIn;
            yb = coefp*xp;
            yp = yb;
            while (x--)
            {
                f32 yc = a0*(f32)(*l_pdImgIn) + p_dA1*(s32)xp - b1*yp - b2*yb;
                *l_pus_J++ = (s32)(yc+0.5);
                xp = *l_pdImgIn++;//xc;
                yb = yp;
                yp = yc;
            }
        }
    }
    // reverse pass
    // ensures response is symmetrical
    {
        s16 x;
        s16 y;
        s32 xn;
        s32 xa;
        f32 flYn;
        f32 ya;
        l_pdImgIn = p_pdImgIn+(s32)w-1+(s32)(h-1)*(s32)w;
        l_pus_J = p_pdImgOut+(s32)w-1+(s32)(h-1)*(s32)w;

        for (y=h-1;y>=0;y--)
        {
            xn = xa = *l_pdImgIn;
            flYn = coefn*(f32)xn;
            ya = flYn;
            x = w;
            while (x--)
            {
                f32 yc = a2*(f32)xn + a3*(f32)xa - b1*flYn - b2*ya;
                xa = xn;
                xn = *l_pdImgIn--;//xc;
                ya = flYn;
                flYn = yc;
                *l_pus_J = *l_pus_J + (s32)(yc+0.5);
                l_pus_J--;
            }
        }
    }
}

void c_DataFusion::Gaussian_Recursive_y(s32 *p_pdImgOut,const s32 *p_pdImgIn, s16 w, s16 h,
            const f32 a0, const f32 p_dA1, const f32 a2, const f32 a3,
            const f32 b1, const f32 b2, const f32 coefp, const f32 coefn)
{
    const s32 *l_pdImgIn ;
    s32 *l_pus_J;
    {
        s16 x;
        s16 y;
        f32 yc;
        f32 yp;  // previous output
        f32 yb;  // previous output by 2
        s32 xp;  // previous input
        l_pdImgIn =  p_pdImgIn;
        l_pus_J =  p_pdImgOut;

        for (x=0;x<w;x++)
        {
            l_pus_J = p_pdImgOut+ x;
            l_pdImgIn  = p_pdImgIn+ x;
            xp = *l_pdImgIn;
            yb = coefp*xp;
            yp = yb;
            for (y=0;y<h;y++)
            {
                yc = a0*(f32)(*l_pdImgIn) + p_dA1*(f32)xp - b1*yp - b2*yb;
                *l_pus_J = (s32)(yc+0.5);
                l_pus_J+= w;
                xp = (*l_pdImgIn);
                l_pdImgIn += w;
                yb = yp;
                yp = yc;
            }
        }
    }
    // reverse pass
    // ensures response is symmetrical
    {
        s16 x;
        s16 y;
        s32 xn;
        s32 xa;
        f32 flYn;
        f32 ya;
        f32 yc;

        for (x=0;x<w;x++)
        {
            l_pus_J = p_pdImgOut+ (s32)x + (s32)(h-1)*(s32)w ;
            l_pdImgIn  = p_pdImgIn + (s32)x + (s32)(h-1)*(s32)w;
            xn = xa = *l_pdImgIn;
            flYn = coefn*xn;
            ya = flYn;
            for (y=h-1;y>=0;y--)
            {
                yc = a2*(f32)xn + a3*(f32)xa - b1*flYn - b2*ya;
                xa = xn;
                xn = *l_pdImgIn;
                ya = flYn;
                flYn = yc;
                //*l_pus_J = *l_pus_J + yc;
                *l_pus_J += (s32)(yc+0.5);
                l_pus_J-= w;
                l_pdImgIn -= w;

            }
        }
    }
}

/* Set sigma S
 */
void c_DataFusion::SetSigmaS(float p_fVal)
{    
    m_DataFusion.m_Set.m_fSigma_spatial = p_fVal;
    // Initialize the Gaussian filter
    Gaussian_Recursive_Order0_Init(m_DataFusion.m_Set.m_fSigma_spatial/(f32)(m_DataFusion.m_Set.m_sScaleFactor), m_DataFusion.m_h, m_DataFusion.m_w, &m_DataFusion.m_StGaussian);
}

/* Set sigma R
 */
void c_DataFusion::SetSigmaR(float p_fVal)
{
    m_DataFusion.m_Set.m_fSigma_range = p_fVal;
    // Gaussian LUT for the fusion filters
    FillColorWeightedTable(m_DataFusion.m_ucaTable, m_DataFusion.m_Set.m_fSigma_range, MLF_APP_2D_RANGE);
    // Gaussian LUT for the fusion filters (BF)
    FillColorWeightedTable(m_DataFusion.m_ucaTableBF, m_DataFusion.m_Set.m_fSigma_range, MLF_APP_INTENSITY_RANGE);
}

/* Set sigma Q
 */
void c_DataFusion::SetSigmaQ(float p_fVal)
{
    m_DataFusion.m_Set.m_fSigma_Q = p_fVal;
}

/* Compute credibility map Q
 */
void c_DataFusion::Compute_Q_D()
{
    cv::Mat l_cvD(cv::Size(m_DataFusion.m_w_original, m_DataFusion.m_h_original), CV_16UC1, (unsigned char*)m_DataFusion.m_pt_D, cv::Mat::AUTO_STEP);
    cv::Mat l_cvQ(m_DataFusion.m_w_original, m_DataFusion.m_h_original, CV_8UC1);
    int l_iScale = 1;
    int l_iDelta = 0;
    int l_iDDepth = CV_16S;
    double l_dMinVal, l_dMaxVal;

    // Generate grad_x and grad_y
    cv::Mat l_cvGrad_x, l_cvGrad_y;
    //cv::Mat l_cvAbsGrad_x, l_cvAbsGrad_y;
    cv::Mat l_cvGradient(cv::Size(m_DataFusion.m_w_original, m_DataFusion.m_h_original), CV_16UC1);
    cv::Mat l_cvAbsGradient;

    cv::minMaxLoc(l_cvD, &l_dMinVal, &l_dMaxVal); // Find minimum and maximum depth values
    //std::cout << "Q_D, l_cvD min: " << l_dMinVal << " max: " << l_dMaxVal << std::endl;

    // Gradient X
    //cv::Scharr(l_cvD, l_cvGrad_x, l_iDDepth, 1, 0, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    cv::Sobel(l_cvD, l_cvGrad_x, l_iDDepth, 1, 0, 3, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    l_cvGrad_x = cv::abs(l_cvGrad_x);
    cv::minMaxLoc(l_cvGrad_x, &l_dMinVal, &l_dMaxVal); // Find minimum and maximum depth values
    //l_cvAbsGrad_x = l_cvGrad_x.clone();
    //l_cvAbsGrad_x.convertTo(l_cvAbsGrad_x, CV_8U, 255.0/(l_dMaxVal-l_dMinVal), -l_dMinVal*255.0/(l_dMaxVal-l_dMinVal));
    //std::cout << "Q_D, l_cvGrad_x min: " << l_dMinVal << " max: " << l_dMaxVal << std::endl;

    // Gradient Y
    //cv::Scharr(l_cvD, l_cvGrad_y, l_iDDepth, 0, 1, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    cv::Sobel(l_cvD, l_cvGrad_y, l_iDDepth, 0, 1, 3, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    l_cvGrad_y = cv::abs(l_cvGrad_y);
    cv::minMaxLoc(l_cvGrad_y, &l_dMinVal, &l_dMaxVal); // Find minimum and maximum depth values
    //l_cvAbsGrad_y = l_cvGrad_y.clone();
    //l_cvAbsGrad_y.convertTo(l_cvAbsGrad_y, CV_8U, 255.0/(l_dMaxVal-l_dMinVal), -l_dMinVal*255.0/(l_dMaxVal-l_dMinVal));
    // Total Gradient, approx since it should be G = sqrt(G_x^2 + G_y^2)
    //cv::addWeighted(l_cvGrad_x, 0.5, l_cvGrad_y, 0.5, 0, l_cvGradient);
    long l_lSize = l_cvGradient.rows*l_cvGradient.cols;
    u16* l_pusGradX = (u16*)l_cvGrad_x.data;
    u16* l_pusGradY = (u16*)l_cvGrad_y.data;
    u16* l_pusGradient = (u16*)l_cvGradient.data;
    while (l_lSize--)
    {
        // Note: the gradient might be sqrt(l_pusGradX^2 + l_pusGradY^2) but I remove the sqrt() because after I've to weight it by the exponential exp(-grad^2).
        // Thus, no need to make the power of 2 when weighting...
        *l_pusGradient = (u16)(((*l_pusGradX)*(*l_pusGradX)+(*l_pusGradY)*(*l_pusGradY))/16); // I divide by 16 (4^2, since I've previsously removed the sqrt) to normalize the Sobel kernel (Grad_X, Grad_Y)
        l_pusGradient++;
        l_pusGradX++;
        l_pusGradY++;
    }
    // Dilating by downsampling factor
    cv::Mat l_cvGradientDilated;
    s16 l_sScaleFactor = 1<<m_DataFusion.m_nr_shift; // Scale factor 1*2^3
    cv::Mat l_cvElement = cv::getStructuringElement(0, cv::Size(11, 11), cv::Point(5, 5));
    cv::dilate(l_cvGradient, l_cvGradientDilated, l_cvElement);
    u16* l_pusGradientDilated = (u16*)l_cvGradientDilated.data;
    u08* l_pucQ = (u08*)l_cvQ.data;
    u16  l_sAux;
    l_lSize = l_cvGradientDilated.rows*l_cvGradientDilated.cols;
    while (l_lSize--)
    {
        l_sAux = (u16)(exp(-(*l_pusGradientDilated)/(2*m_DataFusion.m_Set.m_fSigma_Q*m_DataFusion.m_Set.m_fSigma_Q))*255); // Weighting of the credibility map
        if (l_sAux > 255)
            *l_pucQ = 255;
        else
            *l_pucQ = l_sAux;
        l_pusGradientDilated++;
        l_pucQ++;
    }
    // Weight the gradient
    //cv::Mat l_cvGaussianKernel = cv::getGaussianKernel(cv::Size(9, 9), 50, CV_32F);

    //std::cout << "Kernel: " << std::endl << l_cvGaussianKernel << std::endl;
    //cv::filter2D(l_cvDilated, l_cvDilated, l_iDDepth, l_cvGaussianKernel, cv::Point(-1, -1), l_iDelta, cv::BORDER_DEFAULT);

    // Gaussian Blur
    //cv::GaussianBlur(l_cvGradient, l_cvQ, cv::Size(9, 9), m_DataFusion.m_Set.m_fSigma_Q, m_DataFusion.m_Set.m_fSigma_Q, cv::BORDER_DEFAULT);
    //cv::minMaxLoc(l_cvQ, &l_dMinVal, &l_dMaxVal); // Find minimum and maximum depth values
    //l_cvQ.convertTo(l_cvQ, CV_8U, 255.0/(l_dMaxVal-l_dMinVal), -l_dMinVal*255.0/(l_dMaxVal-l_dMinVal));
/*
    u16* l_pusD = (u16*)l_cvD.data;
    l_pusGradX = (u16*)l_cvGrad_x.data;
    l_pusGradY = (u16*)l_cvGrad_y.data;
    l_pusGradient = (u16*)l_cvGradient.data;
    for (short j=0; j<l_cvQ.rows; j++)
    {
        for (short i=0; i<l_cvQ.cols; i++)
        {
            if ((i==64) && (j==48)) // Theoretically, there is no gradient
            {
                std::cout << "D(65,49): " << *l_pusD <<std::endl;
                std::cout << "Grad_X(65,49): " << *l_pusGradX << std::endl;
                std::cout << "Grad_Y(65,49): " << *l_pusGradY << std::endl;
                std::cout << "Gradient(65,49): " << *l_pusGradient << std::endl;
            }
            if ((i==284) && (j==83))
            {
                std::cout << "D(285,84): " << *l_pusD << ", D(286,84): " << *(l_pusD+1) << std::endl;
                std::cout << "Grad_X(285,84): " << *l_pusGradX << ", Grad_X(286,84): " << *(l_pusGradX+1) << std::endl;
                std::cout << "Grad_Y(285,84): " << *l_pusGradY << ", Grad_Y(286,84): " << *(l_pusGradY+1) << std::endl;
                std::cout << "Gradient(285,84): " << *l_pusGradient << ", Gradient(286,84): " << *(l_pusGradient+1) << std::endl;
            }
            *l_pucQ++ = 255 - *l_pucQ;
            l_pusD++;
            l_pusGradX++;
            l_pusGradY++;
            l_pusGradient++;
        }
    }*/
    memcpy(m_DataFusion.m_pt_Q, (u08*)l_cvQ.data, D_MEM_u08);
/*
    cv::namedWindow("SobelX", CV_WINDOW_AUTOSIZE);
    cv::imshow("SobelX", l_cvAbsGrad_x);

    cv::namedWindow("SobelY", CV_WINDOW_AUTOSIZE);
    cv::imshow("SobelY", l_cvAbsGrad_y);

    cv::namedWindow("Gradient", CV_WINDOW_AUTOSIZE);
    cv::imshow("Gradient", l_cvAbsGradient);

    cv::namedWindow("GaussianBlur", CV_WINDOW_AUTOSIZE);
    cv::imshow("GaussianBlur", l_cvQ);

    cv::waitKey(0);*/
}

/* Compute credibility map Q for the RGB-D filter
 */
void c_DataFusion::Compute_Q_D_RGBD_Filter()
{
    cv::Mat l_cvD(cv::Size(m_DataFusion.m_w_original, m_DataFusion.m_h_original), CV_16UC1, (unsigned char*)m_DataFusion.m_pt_D, cv::Mat::AUTO_STEP);
    cv::Mat l_cvQ(m_DataFusion.m_w_original, m_DataFusion.m_h_original, CV_8UC1);
    int l_iScale = 1;
    int l_iDelta = 0;
    int l_iDDepth = CV_16S;
    double l_dMinVal, l_dMaxVal;

    // Generate grad_x and grad_y
    cv::Mat l_cvGrad_x, l_cvGrad_y;
    cv::Mat l_cvGradient(cv::Size(m_DataFusion.m_w_original, m_DataFusion.m_h_original), CV_16UC1);
    cv::minMaxLoc(l_cvD, &l_dMinVal, &l_dMaxVal); // Find minimum and maximum depth values

    // Gradient X
    //cv::Scharr(l_cvD, l_cvGrad_x, l_iDDepth, 1, 0, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    cv::Sobel(l_cvD, l_cvGrad_x, l_iDDepth, 1, 0, 3, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    l_cvGrad_x = cv::abs(l_cvGrad_x);
    cv::minMaxLoc(l_cvGrad_x, &l_dMinVal, &l_dMaxVal); // Find minimum and maximum depth values

    // Gradient Y
    //cv::Scharr(l_cvD, l_cvGrad_y, l_iDDepth, 0, 1, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    cv::Sobel(l_cvD, l_cvGrad_y, l_iDDepth, 0, 1, 3, l_iScale, l_iDelta, cv::BORDER_DEFAULT);
    l_cvGrad_y = cv::abs(l_cvGrad_y);
    cv::minMaxLoc(l_cvGrad_y, &l_dMinVal, &l_dMaxVal); // Find minimum and maximum depth values

    long l_lSize = l_cvGradient.rows*l_cvGradient.cols;
    u16* l_pusGradX = (u16*)l_cvGrad_x.data;
    u16* l_pusGradY = (u16*)l_cvGrad_y.data;
    u16* l_pusGradient = (u16*)l_cvGradient.data;
    while (l_lSize--)
    {
        // Note: the gradient might be sqrt(l_pusGradX^2 + l_pusGradY^2) but I remove the sqrt() because after I've to weight it by the exponential exp(-grad^2).
        // Thus, no need to make the power of 2 when weighting...
        *l_pusGradient = (u16)(((*l_pusGradX)*(*l_pusGradX)+(*l_pusGradY)*(*l_pusGradY))/16); // I divide by 16 (4^2, since I've previsously removed the sqrt) to normalize the Sobel kernel (Grad_X, Grad_Y)
        l_pusGradient++;
        l_pusGradX++;
        l_pusGradY++;
    }
    // Dilating by 3x3, since this Credibility Map does not consider low-res depth maps as input
    cv::Mat l_cvGradientDilated;
    s16 l_sScaleFactor = 1<<1; // Scale factor 1*2^1
    cv::Mat l_cvElement = cv::getStructuringElement(0, cv::Size(l_sScaleFactor+1, l_sScaleFactor+1), cv::Point(l_sScaleFactor/2, l_sScaleFactor/2));
    cv::dilate(l_cvGradient, l_cvGradientDilated, l_cvElement);
    u16* l_pusGradientDilated = (u16*)l_cvGradientDilated.data;
    u08* l_pucQ = (u08*)l_cvQ.data;
    u16* l_pusD = (u16*)l_cvD.data;
    u16  l_sAux;
    l_lSize = l_cvGradientDilated.rows*l_cvGradientDilated.cols;
    while (l_lSize--)
    {
        l_sAux = (u16)(exp(-(*l_pusGradientDilated)/(2*m_DataFusion.m_Set.m_fSigma_Q*m_DataFusion.m_Set.m_fSigma_Q))*255); // Weighting of the credibility map
        if (l_sAux > 255)
            *l_pucQ = 255;
        else
            *l_pucQ = l_sAux;
        if (*l_pusD == MLF_APP_MAX_DISTANCE) //FGa, bg pixels to non reliable in order to be filled by right depth values
            *l_pucQ = 0;
        l_pusGradientDilated++;
        l_pucQ++;
        l_pusD++;
    }

    memcpy(m_DataFusion.m_pt_Q, (u08*)l_cvQ.data, D_MEM_u08);
}

/* Compute the blending mask (beta) that combines J_2 with J_3 or D
 */
void c_DataFusion::ComputeBlendingMask()
{
    cv::Mat l_cv_Q_I_R(D_Fusion_2D_3D_2DH, D_Fusion_2D_3D_2DW, CV_8UC1, m_DataFusion.m_pt_I_R, cv::Mat::AUTO_STEP);
    cv::Mat l_cv_Q_I_G(D_Fusion_2D_3D_2DH, D_Fusion_2D_3D_2DW, CV_8UC1, m_DataFusion.m_pt_I_G, cv::Mat::AUTO_STEP);
    cv::Mat l_cv_Q_I_B(D_Fusion_2D_3D_2DH, D_Fusion_2D_3D_2DW, CV_8UC1, m_DataFusion.m_pt_I_B, cv::Mat::AUTO_STEP);
    l_cv_Q_I_R = l_cv_Q_I_R.clone();
    l_cv_Q_I_G = l_cv_Q_I_G.clone();
    l_cv_Q_I_B = l_cv_Q_I_B.clone();
    // Compute the gradient
    cv::morphologyEx(l_cv_Q_I_R, l_cv_Q_I_R, cv::MORPH_GRADIENT, cv::Mat());
    cv::morphologyEx(l_cv_Q_I_G, l_cv_Q_I_G, cv::MORPH_GRADIENT, cv::Mat());
    cv::morphologyEx(l_cv_Q_I_B, l_cv_Q_I_B, cv::MORPH_GRADIENT, cv::Mat());
    // Dilate by 3 px
    cv::Mat l_cvElement(3, 3, CV_8U, cv::Scalar(1));
    cv::dilate(l_cv_Q_I_R, l_cv_Q_I_R, l_cvElement);
    cv::dilate(l_cv_Q_I_G, l_cv_Q_I_G, l_cvElement);
    cv::dilate(l_cv_Q_I_B, l_cv_Q_I_B, l_cvElement);
    // Weighting
    u08* l_pucQ_R = l_cv_Q_I_R.data;
    u08* l_pucQ_G = l_cv_Q_I_G.data;
    u08* l_pucQ_B = l_cv_Q_I_B.data;
    f32 l_f32SigmaSq = (f32)((s32)2*(s32)5*(s32)5); // FGa, todo, the 5 has to be the sigma range to weight the credibility map and may be a parameter
    f32 l_fintVal;
    s32 i = D_Fusion_2D_3D_2DSIZE;

    while(i--)
    {
        /* CM = exp(-x2/signma2) where  x = sqrt(l_l1)/4*/
        // Red Channel
        l_fintVal = exp(-((f32)(*l_pucQ_R)/l_f32SigmaSq));
        l_fintVal *= 255;
        if(l_fintVal > 255)
            *l_pucQ_R++ = 255;
        else if(l_fintVal < 1)  // avoid weight of 0
            *l_pucQ_R++ = (u08)0;
        else
            *l_pucQ_R++ = (u08)l_fintVal;
        // Green Channel
        l_fintVal = exp(-((f32)(*l_pucQ_G)/l_f32SigmaSq));
        l_fintVal *= 255;
        if(l_fintVal > 255)
            *l_pucQ_G++ = 255;
        else if(l_fintVal < 1)  // avoid weight of 0
            *l_pucQ_G++ = (u08)0;
        else
            *l_pucQ_G++ = (u08)l_fintVal;
        // Blue Channel
        l_fintVal = exp(-((f32)(*l_pucQ_B)/l_f32SigmaSq));
        l_fintVal *= 255;
        if(l_fintVal > 255)
            *l_pucQ_B++ = 255;
        else if(l_fintVal < 1)  // avoid weight of 0
            *l_pucQ_B++ = (u08)0;
        else
            *l_pucQ_B++ = (u08)l_fintVal;
    }
    // Compute the Beta function (HR) that combines the interpolated J_2 and D
    u08* l_pt_BetaValue = m_DataFusion.m_pt_BetaValue;
    u08* l_pt_BetaCannel = m_DataFusion.m_pt_BetaChannel;
    u08* l_puc_Q_D = m_DataFusion.m_pt_Q;
    u08 l_ucQ_I;
    u08 l_ucQ_D;
    u16 l_usVal;
    l_pucQ_R = l_cv_Q_I_R.data;
    l_pucQ_G = l_cv_Q_I_G.data;
    l_pucQ_B = l_cv_Q_I_B.data;
    long l_lCounter = D_Fusion_2D_3D_2DSIZE;

    while (l_lCounter--)
    {
        if ((*l_pucQ_R <= *l_pucQ_G) && (*l_pucQ_R <= *l_pucQ_B)) // R-Channel
        {
            l_ucQ_I = *l_pucQ_R;
            *l_pt_BetaCannel = D_CHANNEL_R;
        }
        else if ((*l_pucQ_G <= *l_pucQ_R) && (*l_pucQ_G <= *l_pucQ_B)) // G-Channel
        {
            l_ucQ_I = *l_pucQ_G;
            *l_pt_BetaCannel = D_CHANNEL_G;
        }
        else // B-Channel
        {
            l_ucQ_I = *l_pucQ_B;
            *l_pt_BetaCannel = D_CHANNEL_B;
        }
        l_ucQ_D = *l_puc_Q_D;
        l_usVal = l_ucQ_D + (l_ucQ_I*l_ucQ_D)/255 - ((l_ucQ_I*l_ucQ_D)/255)*l_ucQ_D/255;
        *l_pt_BetaValue = (u08)(l_usVal);
        l_pucQ_R++;
        l_pucQ_G++;
        l_pucQ_B++;
        l_puc_Q_D++;
        l_pt_BetaValue++;
        l_pt_BetaCannel++;
    }
}

/* Downsample grayscale image
 */
void c_DataFusion::Downsample(const u08 *p_pcImgIn, u08 *p_pcImgOut, s16 h, s16 w, s16 scale_exp)
{
    const s32 l_lShift = 1<<scale_exp;
    const s32 l_lShiftIn = l_lShift*w;
    const s16 l_sHo=(h>>scale_exp);
    const s16 l_sWo=(w>>scale_exp);
    const u08 *l_pcImgIn;
    s16 y,x;
    const u08 *l_pcStartIn = p_pcImgIn; // + ((l_lShift/2)*w) + (l_lShift/2); //FGa, to consider the mid pixel of the downsampling window

    for (y=0;y<l_sHo;y++)
    {
        l_pcImgIn = l_pcStartIn + (l_lShiftIn) *  (y);
        x = l_sWo;
        //printf("\nMLF_down_sample_1,  y: %hd ",y);
        while (x--)
        {
            *p_pcImgOut++ =*l_pcImgIn;
            l_pcImgIn +=  l_lShift;
        }
    }
}

/* Downsample depth map
 */
void c_DataFusion::Downsample(const u16 *p_pcImgIn, u16 *p_pcImgOut, s16 h, s16 w, s16 scale_exp)
{
    const s32 l_lShift = 1<<scale_exp;
    const s32 l_lShiftIn = l_lShift*w;
    const s16 l_sHo=(h>>scale_exp);
    const s16 l_sWo=(w>>scale_exp);
    const u16 *l_pcImgIn;
    s16 y,x;
    const u16 *l_pcStartIn = p_pcImgIn;// + ((l_lShift/2)*w) + (l_lShift/2); //FGa, to consider the mid pixel of the downsampling window

    for (y=0;y<l_sHo;y++)
    {
        l_pcImgIn = l_pcStartIn + (l_lShiftIn) *  (y);
        x = l_sWo;
        //printf("\nMLF_down_sample_1,  y: %hd ",y);
        while (x--)
        {            
            *p_pcImgOut++ = *l_pcImgIn;
            l_pcImgIn +=  l_lShift;
        }
    }
}

/* Upsample depth map
 */
void c_DataFusion::Upsample(const u16 *p_psImgIn, u16 *p_psImgOut, s16 h, s16 w, s16 scale_exp)
{
    const s16 l_sHo=(h<<scale_exp); // 480
    const s16 l_sWo=(w<<scale_exp); // 640
    const u16* l_ptImgIn;
    u16* l_ptImgOut;

    for (s16 y=0;y<l_sHo;y++)
    {
        for (s16 x=0;x<l_sWo;x++)
        {
            l_ptImgIn = p_psImgIn + ((s16)(y>>scale_exp))*w +(s16)(x>>scale_exp);
            l_ptImgOut = p_psImgOut + y*l_sWo + x;
            *l_ptImgOut = *l_ptImgIn;
        }
    }
}

/* Get min max values (u08 values)
 */
void c_DataFusion::GetMinMaxVal(u08* p_pucMin, u08* p_pucMax, const u08 *p_pucImg, s32 p_lLen)
{
    *p_pucMin = *p_pucImg;
    *p_pucMax = *p_pucImg++;
    p_lLen--;

    while (p_lLen--)
    {
        if (*p_pucMin > *p_pucImg)
        {
            *p_pucMin = *p_pucImg;
        }
        else
        {
            if (*p_pucMax < *p_pucImg)
            {
                *p_pucMax = *p_pucImg;
            }
        }
        p_pucImg++;
    }
}

/* Get min max values (u16 values)
 */
void c_DataFusion::GetMinMaxVal(u16* p_pusMin, u16* p_pusMax, const u16 *p_pusImg, s32 p_lLen)
{
    *p_pusMin = *p_pusImg;
    *p_pusMax = *p_pusImg++;
    p_lLen--;

    while (p_lLen--)
    {
        if (*p_pusMin > *p_pusImg)
        {
            *p_pusMin = *p_pusImg;
        }
        else
        {
            if ((*p_pusMax < *p_pusImg) && (*p_pusImg != MLF_APP_MAX_DISTANCE))
            {
                *p_pusMax = *p_pusImg;
            }
        }
        p_pusImg++;
    }
}

/* Get min max values
 */
void c_DataFusion::GetMinMaxVal_RGB(u08* p_pucMin, u08* p_pucMax, const u08 *p_pucImg_R, const u08 *p_pucImg_G, const u08 *p_pucImg_B, s32 p_lLen)
{
    *p_pucMin = *p_pucImg_R;
    *p_pucMax = *p_pucImg_R++;
    s32 l_lLen = p_lLen--;
    // Red channel
    while (l_lLen--)
    {
        if (*p_pucMin > *p_pucImg_R)
        {
            *p_pucMin = *p_pucImg_R;
        }
        else
        {
            if (*p_pucMax < *p_pucImg_R)
            {
                *p_pucMax = *p_pucImg_R;
            }
        }
        p_pucImg_R++;
    }
    // Green channel
    l_lLen = p_lLen--;
    while (l_lLen--)
    {
        if (*p_pucMin > *p_pucImg_G)
        {
            *p_pucMin = *p_pucImg_G;
        }
        else
        {
            if (*p_pucMax < *p_pucImg_G)
            {
                *p_pucMax = *p_pucImg_G;
            }
        }
        p_pucImg_G++;
    }
    // Blue channel
    l_lLen = p_lLen--;
    while (l_lLen--)
    {
        if (*p_pucMin > *p_pucImg_B)
        {
            *p_pucMin = *p_pucImg_B;
        }
        else
        {
            if(*p_pucMax < *p_pucImg_B)
            {
                *p_pucMax = *p_pucImg_B;
            }
        }
        p_pucImg_B++;
    }
}

/* Convert to Grayscale
 */
void c_DataFusion::ConvertToGrayscale(const unsigned char* p_ptImgIn, long p_lImgSize, unsigned char* p_ptImgOutC)
{
    double l_dAux=0;

    while (p_lImgSize--)
    {
        l_dAux = 0.299*(*p_ptImgIn++);
        l_dAux += 0.587*(*p_ptImgIn++);
        l_dAux += 0.114*(*p_ptImgIn++);
        *p_ptImgOutC++ = (unsigned char)(l_dAux);
        //*p_ptImgOutC++ = (unsigned char)(l_lAux /3); // Approx.
    }
}

/* Returns the guidance image I_RGB
 */
u08* c_DataFusion::GetGuidanceImageRGB()
{
    return m_DataFusion.m_pt_I_RGB;
}

/* Returns the guidance image I
 */
u08* c_DataFusion::GetGuidanceImage()
{
    return m_DataFusion.m_pt_I_GRAY;
}

/* Returns the downsampled guidance image I
 */
u08* c_DataFusion::GetGuidanceImage_ds()
{
    return m_DataFusion.m_pt_I_GRAY_DS;
}

/* Returns the input depth map D
 */
u16* c_DataFusion::GetDepthMap()
{
   return m_DataFusion.m_pt_D;
}

/* Returns the downsampled input depth map D
 */
u16* c_DataFusion::GetDepthMap_ds()
{
    return m_DataFusion.m_pt_D_DS;
}

/* Returns the credibility map Q
 */
u08* c_DataFusion::GetCredibilityMap()
{
   return m_DataFusion.m_pt_Q;
}

/* Returns the downsampled credibility map Q
 */
u08* c_DataFusion::GetCredibilityMap_ds()
{
   return m_DataFusion.m_pt_Q_DS;
}

/* Returns the blending mask B
 */
u08* c_DataFusion::GetBlendingMask()
{
    return m_DataFusion.m_pt_BetaChannel;
}

/* Returns the enhanced depth data (after filtering)
 */
u16* c_DataFusion::GetEnhancedDepthData()
{
    return m_DataFusion.m_pt_J;
}

/* Data processing and internal strcuture update
 */
void c_DataFusion::DataProcessing(u16* p_usDepthData, u08* p_ucRGBData, s16 p_usFIlterType)
{
    long l_lCounter = D_Fusion_2D_3D_2DSIZE;
    const u08* l_puc_I = p_ucRGBData;
    const u16* l_pus_D = p_usDepthData;
    u16* l_pt_D = m_DataFusion.m_pt_D;
    u08* l_pt_I_RGB = m_DataFusion.m_pt_I_RGB;
    u08* l_pt_I_R = m_DataFusion.m_pt_I_R;
    u08* l_pt_I_G = m_DataFusion.m_pt_I_G;
    u08* l_pt_I_B = m_DataFusion.m_pt_I_B;

    switch (p_usFIlterType)
    {
    case BF_FILTER:
        // Copy the input data to the internal data structure        
        while (l_lCounter--)
        {
            if (*l_pus_D == 0) // Occluded pixels are setted to maximum distance
                *l_pt_D = MLF_APP_MAX_DISTANCE;
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Upsampling of the depth map D_DS
        Upsample(m_DataFusion.m_pt_D_DS, m_DataFusion.m_pt_D, m_DataFusion.m_h, m_DataFusion.m_w, m_DataFusion.m_nr_shift);
        break;
    case JBU_FILTER:
        // Copy the input data to the internal data structure
        while (l_lCounter--)
        {           
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            if (*l_pus_D == 0) // Occluded pixels are setted to maximum distance
                *l_pt_D = MLF_APP_MAX_DISTANCE;
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Upsampling of the depth map D_DS
        //Upsample(m_DataFusion.m_pt_D_DS, m_DataFusion.m_pt_D, m_DataFusion.m_h, m_DataFusion.m_w, m_DataFusion.m_nr_shift);
        // Convert the guidance image from RGB to grayscale
        ConvertToGrayscale(m_DataFusion.m_pt_I_RGB, D_Fusion_2D_3D_2DSIZE, m_DataFusion.m_pt_I_GRAY);
        // Downsampling of the guidance image I (grayscale)
        Downsample(m_DataFusion.m_pt_I_GRAY, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        break;
    case PWAS_FILTER:
        // Copy the input data to the internal data structure
        while (l_lCounter--)
        {
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            if (*l_pus_D == 0) // Occluded pixels are setted to maximum distance
                *l_pt_D = MLF_APP_MAX_DISTANCE;
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }        
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Upsampling of the depth map D_DS
        //Upsample(m_DataFusion.m_pt_D_DS, m_DataFusion.m_pt_D, m_DataFusion.m_h, m_DataFusion.m_w, m_DataFusion.m_nr_shift);
        // Convert the guidance image from RGB to grayscale
        ConvertToGrayscale(m_DataFusion.m_pt_I_RGB, D_Fusion_2D_3D_2DSIZE, m_DataFusion.m_pt_I_GRAY);
        // Downsampling of the guidance image I (grayscale)
        Downsample(m_DataFusion.m_pt_I_GRAY, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the credibility map Q (occluded pixels might be preserved)
        Compute_Q_D();
        // Downsampling of the credibility map Q
        Downsample(m_DataFusion.m_pt_Q, m_DataFusion.m_pt_Q_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        break;
    case UML_FILTER:
        // Copy the input data to the internal data structure
        while (l_lCounter--)
        {
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            if (*l_pus_D == 0) // Occluded pixels are setted to maximum distance
                *l_pt_D = MLF_APP_MAX_DISTANCE;
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Upsampling of the depth map D_DS
        // Upsample(m_DataFusion.m_pt_D_DS, m_DataFusion.m_pt_D, m_DataFusion.m_h, m_DataFusion.m_w, m_DataFusion.m_nr_shift);        
        // Compute the credibility map Q (occluded pixels might be preserved)
        Compute_Q_D();
        // Downsampling of the credibility map Q
        Downsample(m_DataFusion.m_pt_Q, m_DataFusion.m_pt_Q_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the Bilateral Filter
        // Initialize the Gaussian filter
        Gaussian_Recursive_Order0_Init(5/(f32)(m_DataFusion.m_Set.m_sScaleFactor), m_DataFusion.m_h, m_DataFusion.m_w, &m_DataFusion.m_StGaussian);
        FillColorWeightedTable(m_DataFusion.m_ucaTableBF, 100, MLF_APP_INTENSITY_RANGE);
        BF_Filter();
        memcpy(m_DataFusion.m_pt_D, m_DataFusion.m_pt_J, D_MEM_u16);
        //Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Convert the guidance image from RGB to grayscale
        ConvertToGrayscale(m_DataFusion.m_pt_I_RGB, D_Fusion_2D_3D_2DSIZE, m_DataFusion.m_pt_I_GRAY);
        // Downsampling of the guidance image I (grayscale)
        Downsample(m_DataFusion.m_pt_I_GRAY, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        Gaussian_Recursive_Order0_Init(m_DataFusion.m_Set.m_fSigma_spatial/(f32)(m_DataFusion.m_Set.m_sScaleFactor), m_DataFusion.m_h, m_DataFusion.m_w, &m_DataFusion.m_StGaussian);
        break;
    case PWAS_RGB_FILTER:
        // Copy the input data to the internal data structure
        while (l_lCounter--)
        {
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_R++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_G++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_B++ = *l_puc_I++;
            if (*l_pus_D == 0) // Occluded pixels are setted to maximum distance
                *l_pt_D = MLF_APP_MAX_DISTANCE;
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Upsampling of the depth map D_DS
        Upsample(m_DataFusion.m_pt_D_DS, m_DataFusion.m_pt_D, m_DataFusion.m_h, m_DataFusion.m_w, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_R
        Downsample(m_DataFusion.m_pt_I_R, m_DataFusion.m_pt_I_R_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_G
        Downsample(m_DataFusion.m_pt_I_G, m_DataFusion.m_pt_I_G_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_B
        Downsample(m_DataFusion.m_pt_I_B, m_DataFusion.m_pt_I_B_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Convert the guidance image from RGB to grayscale
        ConvertToGrayscale(m_DataFusion.m_pt_I_RGB, D_Fusion_2D_3D_2DSIZE, m_DataFusion.m_pt_I_GRAY);
        // Downsampling of the guidance image I (grayscale)
        Downsample(m_DataFusion.m_pt_I_GRAY, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the credibility map Q (occluded pixels might be preserved)
        Compute_Q_D();
        // Downsampling of the credibility map Q
        Downsample(m_DataFusion.m_pt_Q, m_DataFusion.m_pt_Q_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the Beta function (HR) that combines the interpolated J_2 and D (GRAY channel not used!)
        ComputeBlendingMask();
        break;
    case UML_RGB_FILTER:
        // Copy the input data to the internal data structure
        while (l_lCounter--)
        {
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_R++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_G++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_B++ = *l_puc_I++;
            if (*l_pus_D == 0) // Occluded pixels are setted to maximum distance
                *l_pt_D = MLF_APP_MAX_DISTANCE;
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Upsampling of the depth map D_DS
        Upsample(m_DataFusion.m_pt_D_DS, m_DataFusion.m_pt_D, m_DataFusion.m_h, m_DataFusion.m_w, m_DataFusion.m_nr_shift);
        // Compute the Bilateral Filter
        BF_Filter();
        memcpy(m_DataFusion.m_pt_D, m_DataFusion.m_pt_J, D_MEM_u16);
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_R
        Downsample(m_DataFusion.m_pt_I_R, m_DataFusion.m_pt_I_R_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_G
        Downsample(m_DataFusion.m_pt_I_G, m_DataFusion.m_pt_I_G_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_B
        Downsample(m_DataFusion.m_pt_I_B, m_DataFusion.m_pt_I_B_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Convert the guidance image from RGB to grayscale
        ConvertToGrayscale(m_DataFusion.m_pt_I_RGB, D_Fusion_2D_3D_2DSIZE, m_DataFusion.m_pt_I_GRAY);
        // Downsampling of the guidance image I (grayscale)
        Downsample(m_DataFusion.m_pt_I_GRAY, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the credibility map Q (occluded pixels might be preserved)
        Compute_Q_D();
        // Downsampling of the credibility map Q
        Downsample(m_DataFusion.m_pt_Q, m_DataFusion.m_pt_Q_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the Beta function (HR) that combines the interpolated J_2 and D (GRAY channel not used!)
        ComputeBlendingMask();
        break;
    case RGBD_FILTER:
        // Copy the input data to the internal data structure
        while (l_lCounter--)
        {
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_R++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_G++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I;
            *l_pt_I_B++ = *l_puc_I++;
            if (*l_pus_D == 0) // Occluded pixels are setted to maximum distance
                *l_pt_D = MLF_APP_MAX_DISTANCE;
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_R
        Downsample(m_DataFusion.m_pt_I_R, m_DataFusion.m_pt_I_R_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_G
        Downsample(m_DataFusion.m_pt_I_G, m_DataFusion.m_pt_I_G_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Downsampling the guidance image I_B
        Downsample(m_DataFusion.m_pt_I_B, m_DataFusion.m_pt_I_B_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Convert the guidance image from RGB to grayscale
        ConvertToGrayscale(m_DataFusion.m_pt_I_RGB, D_Fusion_2D_3D_2DSIZE, m_DataFusion.m_pt_I_GRAY);
        // Downsampling of the guidance image I (grayscale)
        Downsample(m_DataFusion.m_pt_I_GRAY, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the credibility map Q (occluded pixels might be filled)
        Compute_Q_D_RGBD_Filter();
        // Downsampling of the credibility map Q
        Downsample(m_DataFusion.m_pt_Q, m_DataFusion.m_pt_Q_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Compute the Beta function (HR) that combines the interpolated J_2 and D (GRAY channel not used!)
        ComputeBlendingMask();
        break;
    case LITERATURE_FILTERS:
        // Copy the input data to the internal data structure
        while (l_lCounter--)
        {
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            *l_pt_I_RGB++ = *l_puc_I++;
            if (*l_pus_D == MLF_APP_MAX_DISTANCE) // Occluded pixels are setted to maximum distance
                *l_pt_D = 0; // FGa, bg pixels are setted to 0 in order to less influence when dilating to compute the max diff omega
            else
                *l_pt_D = *l_pus_D;
            l_pt_D++;
            l_pus_D++;
        }
        // Downsampling of the depth map D
        Downsample(m_DataFusion.m_pt_D, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        // Upsampling of the depth map D_DS
        Upsample(m_DataFusion.m_pt_D_DS, m_DataFusion.m_pt_D, m_DataFusion.m_h, m_DataFusion.m_w, m_DataFusion.m_nr_shift);
        // Convert the guidance image from RGB to grayscale
        ConvertToGrayscale(m_DataFusion.m_pt_I_RGB, D_Fusion_2D_3D_2DSIZE, m_DataFusion.m_pt_I_GRAY);
        // Downsampling of the guidance image I (grayscale)
        Downsample(m_DataFusion.m_pt_I_GRAY, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h_original, m_DataFusion.m_w_original, m_DataFusion.m_nr_shift);
        break;
    }
}

/*** FILTERS ***/

/* Apply Bilateral Filter (BF Filter)
 */
short c_DataFusion::BF_Filter()
{    
    const u16* l_pt_D_DS;
    const u16* l_pt_I_DS;
    u16 l_tImage_min;
    u16 l_tImage_max;
    u16 l_tImg_Range;
    const u16 l_usQuantizationLevels = D_QUANTIZATION_RANGE_3D-1;

    if (!m_DataFusion.m_sInitOk)
    {
        return -1;
    }
    if ((u16)D_QUANTIZATION_RANGE_3D > MLF_APP_INTENSITY_RANGE)
    {
        return -2;
    }
#ifdef	MLF_ANALYSE_TIME_CONSUMPTION
    std::cout << "Time consumption analysis for the BF filter: " << std::endl;
    const clock_t c_start = std::clock();
    short l_sNbIter = 1;
    for (int ii=0;ii<l_sNbIter;ii++)
    {
#endif 
    // Compute the range [min,max] of the 3D image. This range will be divided into D_QUANTIZATION_RANGE_3D levels
    GetMinMaxVal(&l_tImage_min, &l_tImage_max, m_DataFusion.m_pt_D_DS, m_DataFusion.m_h*m_DataFusion.m_w);
    l_tImg_Range = l_tImage_max-l_tImage_min;

    s32 *l_pt_Jk[2];
    s32 *l_pl_Jk;
    s16 jk_0 = 0;
    s16 jk_1 = 1;
    s32 y,x;
    u08  l_ucShiftUnit = (1<<m_DataFusion.m_nr_shift);  // 1..16
#ifdef	MLF_BackGround_Counter
    s32 *l_pt_Jk_BG[2];
    s32 *l_pl_Jk_BG;
    l_pt_Jk_BG[0] = m_DataFusion.m_pt_Jk_level0_R_BG;
    l_pt_Jk_BG[1] = m_DataFusion.m_pt_Jk_level1_R_BG;
#endif
    l_pt_Jk[0] = m_DataFusion.m_pt_Jk_level0_R;
    l_pt_Jk[1] = m_DataFusion.m_pt_Jk_level1_R;

    // For all levels
    for (s16 i=0;i<D_QUANTIZATION_RANGE_3D;i++)
    {
        if (i==0)
        {
            l_pl_Jk=	l_pt_Jk[jk_0];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_0];
#endif
        }
        else
        {
            l_pl_Jk=	l_pt_Jk[jk_1];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_1];
#endif
        }
        {
            s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;
            s32* l_pl_Jk_FG = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pl_Jk_BG2 = l_pl_Jk_BG;
#endif
            // RGB to Grayscale
            s32 l_lIndexIni;
            if (i== 0)
                l_lIndexIni = l_tImage_min;
            else
            {
                if (i == l_usQuantizationLevels)
                    l_lIndexIni = l_tImage_max;
                else
                    l_lIndexIni = ((s32)l_tImage_min +(i*(s32)(l_tImg_Range))/l_usQuantizationLevels);
            }
            l_pt_D_DS = m_DataFusion.m_pt_D_DS;
            l_pt_I_DS = m_DataFusion.m_pt_D_DS; // Both spatial and range terms apply to the same source of data (depth)

            for (y=0;y<m_DataFusion.m_h;y++)
            {
                for (x=0;x<m_DataFusion.m_w;x++)
                {
                    s32 index=(l_lIndexIni-(s32)(*l_pt_I_DS++));
                    if (index < 0)
                    {
                        index = -index;
                    }
                    //if (index >= MLF_APP_INTENSITY_RANGE) // Updated by BMi
                    //{
                    //    index = MLF_APP_INTENSITY_RANGE-1;
                    //}
#ifdef	MLF_BackGround_Counter
                    if (*l_pt_D_DS >= MLF_APP_MAX_DISTANCE)
                    {
                        *l_pl_Wk++ = 0x0;
                        *l_pl_Jk_FG++ = 0x0;
                        *l_pl_Jk_BG2++ = (s32)(m_DataFusion.m_ucaTableBF[index]*255); // 255 Increase Accuracy
                    }
                    else
                    {
                        *l_pl_Jk_BG2++ = 0;
                        *l_pl_Wk = (s32)(m_DataFusion.m_ucaTableBF[index]*255); // 255 Increase Accuracy
                        *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                        l_pl_Wk++;
                    }
#else
                    *l_pl_Wk = (s32)(m_DataFusion.m_ucaTableBF[index]*255); // 255 Increase Accuracy
                    *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                    l_pl_Wk++;
#endif
                    l_pt_D_DS++;
                }
            }
        }
        {
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#ifdef	MLF_BackGround_Counter
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#endif
        }
        {
            const s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;

            s32* l_pus_J = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pus_J_BG = l_pl_Jk_BG;
#endif
            y = (s32)m_DataFusion.m_h*(s32)m_DataFusion.m_w;
            while (y--)
            {
                {
#ifdef	MLF_BackGround_Counter
                    if (*l_pus_J_BG == 0)
                        *l_pus_J_BG = 0;
                    else if (*l_pl_Wk <= *l_pus_J_BG)
                        *l_pus_J_BG = 255;
                    else
                        *l_pus_J_BG = (255*(*l_pus_J_BG))/(2*(*l_pl_Wk));
                    l_pus_J_BG++;
#endif
                    if (*l_pl_Wk)
                    {
                        s32 l_lWkHalf = *l_pl_Wk/2; // FGa, to avoid numeric problems (in the test was between 999 and 1000)
                        *l_pus_J = ((*l_pus_J+l_lWkHalf))/(*l_pl_Wk);
                    }
                    l_pus_J++;
                    l_pl_Wk++;
                }
            }
        }
        // INTERPOLATION: Check which pixels in the high-resolution 2D image have an intensity value that belongs to the range (levels 0-1) that we've processed
        //                in order to update the final pixels
        if (i>0)
        {
            const s32 l_lXMax = (s32)((m_DataFusion.m_w-1)*l_ucShiftUnit);
            const s32 l_lYMax = (s32)((m_DataFusion.m_h-1)*l_ucShiftUnit)-1;
#ifdef	MLF_BackGround_Counter
            const s32* p_plImgJk0_BG = l_pt_Jk_BG[jk_0];
            const s32* p_plImgJk0_d_BG = p_plImgJk0_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_BG = l_pt_Jk_BG[jk_1];
            const s32* p_plImgJk1_d_BG = p_plImgJk1_BG+m_DataFusion.m_w;
#endif
            const s32* p_plImgJk0 = l_pt_Jk[jk_0];
            const s32* p_plImgJk0_d = p_plImgJk0+m_DataFusion.m_w;
            const s32* p_plImgJk1 = l_pt_Jk[jk_1];
            const s32* p_plImgJk1_d = p_plImgJk1+m_DataFusion.m_w;
            const s16* l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[0][0][0];
            const s32 l_lImgRangeXSihiftUnitSq = (s32)l_tImg_Range* (s32)l_ucShiftUnit*(s32)l_ucShiftUnit;
            const u16 *l_pt_D = m_DataFusion.m_pt_D;
            const u16 l_tKf1Level = (((i-1)*(s16)l_tImg_Range)/(s16)l_usQuantizationLevels)+l_tImage_min;
            const u16 l_tKf2Level = (((i  )*(s16)l_tImg_Range)/(s16)l_usQuantizationLevels)+l_tImage_min;
            s16 l_sY0 = 0;
            s16 l_sdY = 0;
            u16* l_pt_J = m_DataFusion.m_pt_J;

            for (y=0;y<m_DataFusion.m_h_original;y++)
            {
                {
                    s16 l_sdX = 0;
                    u16 l_usVal;

                    for (x=0;x<l_lXMax;x++)
                    {
                        l_usVal = *l_pt_D;
                        if (l_usVal<l_tKf2Level)
                        {
                           if (l_usVal>=l_tKf1Level)
                           {
                                // qx_linear_interpolate_xy
                                s32 l_lRetInterpolate2;
                                s32 l_lRetInterpolate;
                                const s16 alpha= (l_tKf2Level - l_usVal)* l_usQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d_BG+1));
                                l_lRetInterpolate *= alpha;
                                l_psaDxx -= 3; // reset the position
                                l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate2 *= alpha_Comp;
                                l_psaDxx -= 3; // reset the position
                                l_BakGroundWeight =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d+1));
                                    l_lRetInterpolate *= alpha;

                                    l_psaDxx -= 3;// reset the position
                                    l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate2 *= alpha_Comp;

                                    *l_pt_J =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                }
                            }
                        }
                        else
                        {
                            if (i==l_usQuantizationLevels)
                            {
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
                                s32 l_lRetInterpolate ;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2;
                                l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                l_psaDxx -= 3;// reset the position
                                l_BakGroundWeight = l_lRetInterpolate;
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2;
                                    l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                    *l_pt_J =(u16)(l_lRetInterpolate);
                                }
                            }
                        }
                        l_pt_J++;
                        l_pt_D++;
                        {
                            l_sdX++;
                            if (l_sdX == l_ucShiftUnit)
                            {
                                l_sdX= 0;
#ifdef	MLF_BackGround_Counter
                                p_plImgJk0_BG++;
                                p_plImgJk0_d_BG++;
                                p_plImgJk1_BG++;
                                p_plImgJk1_d_BG++;
#endif
                                p_plImgJk0++;
                                p_plImgJk0_d++;
                                p_plImgJk1++;
                                p_plImgJk1_d++;
                            }
                        }
                    }
                    // X Max to end
                    {
                        s32 l_dLIc;
                        s32 l_dLIa;
                        s32 l_dLIb;
                        const s16* l_psaDxx = l_spaDxx_Ylevel;
                        // qx_linear_interpolate_xy
#ifdef	MLF_BackGround_Counter
                        s32 l_dLIc_BG;
                        s32 l_dLIa_BG;
                        s32 l_dLIb_BG;
                        l_dLIa_BG = ((*l_psaDxx)  *(*(p_plImgJk0_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d_BG)));
                        l_dLIb_BG = (  (*l_psaDxx)*(*(p_plImgJk1_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d_BG)));
                        l_dLIc_BG = l_dLIb_BG+(l_ucShiftUnit*l_ucShiftUnit)/2;
#endif
                        l_dLIa = ((*l_psaDxx)  *(*(p_plImgJk0))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d)));
                        l_dLIb = (  (*l_psaDxx)*(*(p_plImgJk1))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d)));
                        l_dLIc = l_dLIb+(l_ucShiftUnit*l_ucShiftUnit)/2;
                        for (x = l_lXMax;x<m_DataFusion.m_w_original;x++)
                        {
                            if (l_usVal<l_tKf2Level)
                            {
                              if (l_usVal>=l_tKf1Level)
                              {
                                const s16 alpha= (l_tKf2Level - l_usVal)* l_usQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);

#ifdef	MLF_BackGround_Counter
                                *l_pt_J =(u16)(((alpha*l_dLIa_BG+(alpha_Comp)* l_dLIb_BG))/(l_lImgRangeXSihiftUnitSq));
                                if (*l_pt_J > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                    *l_pt_J =(u16)(((alpha*l_dLIa+(alpha_Comp)* l_dLIb))/(l_lImgRangeXSihiftUnitSq));
                              }
                            }
                            else
                            {
                                if (i==l_usQuantizationLevels)
                                {
#ifdef	MLF_BackGround_Counter
                                    *l_pt_J =(u16)(l_dLIc_BG/(l_ucShiftUnit*l_ucShiftUnit));
                                    if (*l_pt_J > D_BG_THRESHOLD)
                                        *l_pt_J = MLF_APP_MAX_DISTANCE;
                                    else
#endif
                                        *l_pt_J =(u16)(l_dLIc/(l_ucShiftUnit*l_ucShiftUnit));
                                }
                            }
                            l_pt_J++;
                            l_pt_D++;
                        }
                    }
                }
                if (y < l_lYMax)
                {
                    l_sdY++;
                    if (l_sdY == l_ucShiftUnit)
                    {
                        l_sdY= 0;
                        l_sY0++;
                    }
                    l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[l_sdY][0][0];
                }
#ifdef	MLF_BackGround_Counter
                p_plImgJk0_BG   = l_pt_Jk_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_BG = p_plImgJk0_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_BG   = l_pt_Jk_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_BG = p_plImgJk1_BG + (s32)m_DataFusion.m_w;
#endif

                p_plImgJk0   = l_pt_Jk[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d = p_plImgJk0 + (s32)m_DataFusion.m_w;
                p_plImgJk1   = l_pt_Jk[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d = p_plImgJk1 + (s32)m_DataFusion.m_w;
            }
            jk_1=jk_0;
            jk_0=(jk_0+1)%2;
        }
    }
#ifdef MLF_ANALYSE_TIME_CONSUMPTION
    }
    std::clock_t c_end = std::clock();

    std::cout << "CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC
              << " ms" << std::endl;
    std::cout << "Done..." << std::endl;
#endif

    return 0;
}

/* Apply Joint Bilateral Upsampling Filter (JBU Filter)
 */
short c_DataFusion::JBU_Filter()
{
    const u16* l_pt_D_DS;
    const u08* l_pt_I_DS;
    u08 l_tImage_min;
    u08 l_tImage_max;
    u08 l_tImg_Range;
    const u08 l_ucQuantizationLevels = D_QUANTIZATION_RANGE_2D-1;

    if (!m_DataFusion.m_sInitOk)
    {
        return -1;
    }
    if ((u16)D_QUANTIZATION_RANGE_2D > MLF_APP_INTENSITY_RANGE)
    {
        return -2;
    }
#ifdef	MLF_ANALYSE_TIME_CONSUMPTION
    std::cout << "Time consumption analysis for the JBU filter: " << std::endl;
    const clock_t c_start = std::clock();
    short l_sNbIter = 1;
    for (int ii=0;ii<l_sNbIter;ii++)
    {
#endif
    // Compute the range [min,max] of the 2D image. This range will be divided into D_QUANTIZATION_RANGE_2D levels
    GetMinMaxVal(&l_tImage_min, &l_tImage_max, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h*m_DataFusion.m_w);
    l_tImg_Range = l_tImage_max-l_tImage_min;
    s32 *l_pt_Jk[2];
    s32 *l_pl_Jk;
    s16 jk_0 = 0;
    s16 jk_1 = 1;
    s32 y,x;
    u08  l_ucShiftUnit = (1<<m_DataFusion.m_nr_shift);  // 1..16
#ifdef	MLF_BackGround_Counter
    s32 *l_pt_Jk_BG[2];
    s32 *l_pl_Jk_BG;
    l_pt_Jk_BG[0] = m_DataFusion.m_pt_Jk_level0_R_BG;
    l_pt_Jk_BG[1] = m_DataFusion.m_pt_Jk_level1_R_BG;
#endif
    l_pt_Jk[0] = m_DataFusion.m_pt_Jk_level0_R;
    l_pt_Jk[1] = m_DataFusion.m_pt_Jk_level1_R;

    // For all levels
    for (s16 i=0;i<D_QUANTIZATION_RANGE_2D;i++)
    {
        if (i==0)
        {
            l_pl_Jk=	l_pt_Jk[jk_0];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_0];
#endif
        }
        else
        {
            l_pl_Jk=	l_pt_Jk[jk_1];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_1];
#endif
        }
        {
            s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;
            s32* l_pl_Jk_FG = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pl_Jk_BG2 = l_pl_Jk_BG;
#endif
            // RGB to Grayscale
            s32 l_lIndexIni;
            if (i== 0)
                l_lIndexIni = l_tImage_min;
            else
            {
                if (i == l_ucQuantizationLevels)
                    l_lIndexIni = l_tImage_max;
                else
                    l_lIndexIni = ((s32)l_tImage_min +(i*(s32)(l_tImg_Range))/l_ucQuantizationLevels);
            }
            l_pt_D_DS = m_DataFusion.m_pt_D_DS;
            l_pt_I_DS = m_DataFusion.m_pt_I_GRAY_DS;

            for (y=0;y<m_DataFusion.m_h;y++)
            {
                for (x=0;x<m_DataFusion.m_w;x++)
                {
                    s32 index=(l_lIndexIni-(s32)(*l_pt_I_DS++));
                    if (index < 0)
                    {
                        index = -index;
                    }
#ifdef	MLF_BackGround_Counter
                    if (*l_pt_D_DS >= MLF_APP_MAX_DISTANCE)
                    {
                        *l_pl_Wk++ = 0x0;
                        *l_pl_Jk_FG++ = 0x0;
                        *l_pl_Jk_BG2++ = (s32)(m_DataFusion.m_ucaTable[index]*255); // 255 Increase Accuracy
                    }
                    else
                    {
                        *l_pl_Jk_BG2++ = 0;
                        *l_pl_Wk = (s32)(m_DataFusion.m_ucaTable[index]*255); // 255 Increase Accuracy
                        *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                        l_pl_Wk++;
                    }
#else
                    *l_pl_Wk = (s32)(m_DataFusion.m_ucaTable[index]*255); // 255 Increase Accuracy
                    *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                    l_pl_Wk++;
#endif
                    l_pt_D_DS++;
                }
            }
        }
        {
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#ifdef	MLF_BackGround_Counter
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#endif
        }
        {
            const s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;

            s32* l_pus_J = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pus_J_BG = l_pl_Jk_BG;
#endif
            y = (s32)m_DataFusion.m_h*(s32)m_DataFusion.m_w;
            while (y--)
            {
                {
#ifdef	MLF_BackGround_Counter
                    if (*l_pus_J_BG == 0)
                        *l_pus_J_BG = 0;
                    else if (*l_pl_Wk <= *l_pus_J_BG)
                        *l_pus_J_BG = 255;
                    else
                        *l_pus_J_BG = (255*(*l_pus_J_BG))/(2*(*l_pl_Wk));
                    l_pus_J_BG++;
#endif
                    if (*l_pl_Wk)
                    {
                        s32 l_lWkHalf = *l_pl_Wk/2; // FGa, to avoid numeric problems (in the test was between 999 and 1000)
                        *l_pus_J = ((*l_pus_J+l_lWkHalf))/(*l_pl_Wk);
                    }
                    l_pus_J++;
                    l_pl_Wk++;
                }
            }
        }
        // INTERPOLATION: Check which pixels in the high-resolution 2D image have an intensity value that belongs to the range (levels 0-1) that we've processed
        //                in order to update the final pixels
        if (i>0)
        {
            const s32 l_lXMax = (s32)((m_DataFusion.m_w-1)*l_ucShiftUnit);
            const s32 l_lYMax = (s32)((m_DataFusion.m_h-1)*l_ucShiftUnit)-1;
#ifdef	MLF_BackGround_Counter
            const s32* p_plImgJk0_BG = l_pt_Jk_BG[jk_0];
            const s32* p_plImgJk0_d_BG = p_plImgJk0_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_BG = l_pt_Jk_BG[jk_1];
            const s32* p_plImgJk1_d_BG = p_plImgJk1_BG+m_DataFusion.m_w;
#endif
            const s32* p_plImgJk0 = l_pt_Jk[jk_0];
            const s32* p_plImgJk0_d = p_plImgJk0+m_DataFusion.m_w;
            const s32* p_plImgJk1 = l_pt_Jk[jk_1];
            const s32* p_plImgJk1_d = p_plImgJk1+m_DataFusion.m_w;
            const s16* l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[0][0][0];
            const s32 l_lImgRangeXSihiftUnitSq = (s32)l_tImg_Range* (s32)l_ucShiftUnit*(s32)l_ucShiftUnit;
            const u08 *l_pt_I = m_DataFusion.m_pt_I_GRAY;
            const u16 l_tKf1Level = (((i-1)*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            const u16 l_tKf2Level = (((i  )*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            s16 l_sY0 = 0;
            s16 l_sdY = 0;
            u16* l_pt_J = m_DataFusion.m_pt_J;

            for (y=0;y<m_DataFusion.m_h_original;y++)
            {
                {
                    s16 l_sdX = 0;
                    u08 l_ucVal;

                    for (x=0;x<l_lXMax;x++)
                    {
                        l_ucVal = *l_pt_I;
                        if (l_ucVal<l_tKf2Level)
                        {
                           if (l_ucVal>=l_tKf1Level)
                           {
                                // qx_linear_interpolate_xy
                                s32 l_lRetInterpolate2;
                                s32 l_lRetInterpolate;
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d_BG+1));
                                l_lRetInterpolate *= alpha;
                                l_psaDxx -= 3; // reset the position
                                l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate2 *= alpha_Comp;
                                l_psaDxx -= 3; // reset the position
                                l_BakGroundWeight =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d+1));
                                    l_lRetInterpolate *= alpha;

                                    l_psaDxx -= 3;// reset the position
                                    l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate2 *= alpha_Comp;

                                    *l_pt_J =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                }
                            }
                        }
                        else
                        {
                            if (i==l_ucQuantizationLevels)
                            {
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
                                s32 l_lRetInterpolate ;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2;
                                l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                l_psaDxx -= 3;// reset the position
                                l_BakGroundWeight = l_lRetInterpolate;
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2;
                                    l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                    *l_pt_J =(u16)(l_lRetInterpolate);
                                }
                            }
                        }
                        l_pt_J++;
                        l_pt_I++;
                        {
                            l_sdX++;
                            if (l_sdX == l_ucShiftUnit)
                            {
                                l_sdX= 0;
#ifdef	MLF_BackGround_Counter
                                p_plImgJk0_BG++;
                                p_plImgJk0_d_BG++;
                                p_plImgJk1_BG++;
                                p_plImgJk1_d_BG++;
#endif

                                p_plImgJk0++;
                                p_plImgJk0_d++;
                                p_plImgJk1++;
                                p_plImgJk1_d++;
                            }
                        }
                    }
                    // X Max to end
                    {
                        s32 l_dLIc;
                        s32 l_dLIa;
                        s32 l_dLIb;
                        const s16* l_psaDxx = l_spaDxx_Ylevel;
                        // qx_linear_interpolate_xy
#ifdef	MLF_BackGround_Counter
                        s32 l_dLIc_BG;
                        s32 l_dLIa_BG;
                        s32 l_dLIb_BG;
                        l_dLIa_BG = ((*l_psaDxx)  *(*(p_plImgJk0_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d_BG)));
                        l_dLIb_BG = (  (*l_psaDxx)*(*(p_plImgJk1_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d_BG)));
                        l_dLIc_BG = l_dLIb_BG+(l_ucShiftUnit*l_ucShiftUnit)/2;
#endif
                        l_dLIa = ((*l_psaDxx)  *(*(p_plImgJk0))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d)));
                        l_dLIb = (  (*l_psaDxx)*(*(p_plImgJk1))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d)));
                        l_dLIc = l_dLIb+(l_ucShiftUnit*l_ucShiftUnit)/2;
                        for (x = l_lXMax;x<m_DataFusion.m_w_original;x++)
                        {
                            if (l_ucVal<l_tKf2Level)
                            {
                              if (l_ucVal>=l_tKf1Level)
                              {
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);

#ifdef	MLF_BackGround_Counter
                                *l_pt_J =(u16)(((alpha*l_dLIa_BG+(alpha_Comp)* l_dLIb_BG))/(l_lImgRangeXSihiftUnitSq));
                                if (*l_pt_J > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                    *l_pt_J =(u16)(((alpha*l_dLIa+(alpha_Comp)* l_dLIb))/(l_lImgRangeXSihiftUnitSq));
                              }
                            }
                            else
                            {
                                if (i==l_ucQuantizationLevels)
                                {
#ifdef	MLF_BackGround_Counter
                                    *l_pt_J =(u16)(l_dLIc_BG/(l_ucShiftUnit*l_ucShiftUnit));
                                    if (*l_pt_J > D_BG_THRESHOLD)
                                        *l_pt_J = MLF_APP_MAX_DISTANCE;
                                    else
#endif
                                        *l_pt_J =(u16)(l_dLIc/(l_ucShiftUnit*l_ucShiftUnit));
                                }
                            }
                            l_pt_J++;
                            l_pt_I++;
                        }
                    }
                }
                if (y < l_lYMax)
                {
                    l_sdY++;
                    if (l_sdY == l_ucShiftUnit)
                    {
                        l_sdY= 0;
                        l_sY0++;
                    }
                    l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[l_sdY][0][0];
                }
#ifdef	MLF_BackGround_Counter
                p_plImgJk0_BG   = l_pt_Jk_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_BG = p_plImgJk0_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_BG   = l_pt_Jk_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_BG = p_plImgJk1_BG + (s32)m_DataFusion.m_w;
#endif

                p_plImgJk0   = l_pt_Jk[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d = p_plImgJk0 + (s32)m_DataFusion.m_w;
                p_plImgJk1   = l_pt_Jk[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d = p_plImgJk1 + (s32)m_DataFusion.m_w;
            }
            jk_1=jk_0;
            jk_0=(jk_0+1)%2;
        }
    }
#ifdef MLF_ANALYSE_TIME_CONSUMPTION
    }
    std::clock_t c_end = std::clock();

    std::cout << "CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC
              << " ms" << std::endl;
    std::cout << "Done..." << std::endl;
#endif

    return 0;
}

/* Apply Pixel Weighted Average Strategy Filter (PWAS Filter)
 */
short c_DataFusion::PWAS_Filter()
{
    const u16* l_pt_D_DS;
    const u08* l_pt_I_DS;
    const u08* l_pt_Q_DS;
    u08 l_tImage_min;
    u08 l_tImage_max;
    u08 l_tImg_Range;
    const u08 l_ucQuantizationLevels = D_QUANTIZATION_RANGE_2D-1;

    if (!m_DataFusion.m_sInitOk)
    {
        return -1;
    }
    if ((u16)D_QUANTIZATION_RANGE_2D > MLF_APP_INTENSITY_RANGE)
    {
        return -2;
    }
#ifdef	MLF_ANALYSE_TIME_CONSUMPTION
    std::cout << "Time consumption analysis for the PWAS filter: " << std::endl;
    const clock_t c_start = std::clock();
    short l_sNbIter = 1;
    for (int ii=0;ii<l_sNbIter;ii++)
    {
#endif
    // Compute the range [min,max] of the 2D image. This range will be divided into D_QUANTIZATION_RANGE_2D levels
    GetMinMaxVal(&l_tImage_min, &l_tImage_max, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h*m_DataFusion.m_w);
    l_tImg_Range = l_tImage_max-l_tImage_min;
    s32 *l_pt_Jk[2];
    s32 *l_pl_Jk;
    s16 jk_0 = 0;
    s16 jk_1 = 1;
    s32 y,x;
    u08  l_ucShiftUnit = (1<<m_DataFusion.m_nr_shift);  // 1..16
#ifdef	MLF_BackGround_Counter
    s32 *l_pt_Jk_BG[2];
    s32 *l_pl_Jk_BG;
    l_pt_Jk_BG[0] = m_DataFusion.m_pt_Jk_level0_R_BG;
    l_pt_Jk_BG[1] = m_DataFusion.m_pt_Jk_level1_R_BG;
#endif
    l_pt_Jk[0] = m_DataFusion.m_pt_Jk_level0_R;
    l_pt_Jk[1] = m_DataFusion.m_pt_Jk_level1_R;

    // For all levels
    for (s16 i=0;i<D_QUANTIZATION_RANGE_2D;i++)
    {
        if (i==0)
        {
            l_pl_Jk=	l_pt_Jk[jk_0];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_0];
#endif
        }
        else
        {
            l_pl_Jk=	l_pt_Jk[jk_1];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_1];
#endif
        }
        {
            s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;
            s32* l_pl_Jk_FG = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pl_Jk_BG2 = l_pl_Jk_BG;
#endif
            // RGB to Grayscale
            s32 l_lIndexIni;
            if (i== 0)
                l_lIndexIni = l_tImage_min;
            else
            {
                if (i == l_ucQuantizationLevels)
                    l_lIndexIni = l_tImage_max;
                else
                    l_lIndexIni = ((s32)l_tImage_min +(i*(s32)(l_tImg_Range))/l_ucQuantizationLevels);
            }
            l_pt_D_DS = m_DataFusion.m_pt_D_DS;
            l_pt_I_DS = m_DataFusion.m_pt_I_GRAY_DS;
            l_pt_Q_DS = m_DataFusion.m_pt_Q_DS;

            for (y=0;y<m_DataFusion.m_h;y++)
            {
                for (x=0;x<m_DataFusion.m_w;x++)
                {
                    s32 index=(l_lIndexIni-(s32)(*l_pt_I_DS++));
                    if (index < 0)
                    {
                        index = -index;
                    }
#ifdef	MLF_BackGround_Counter
                    if (*l_pt_D_DS >= MLF_APP_MAX_DISTANCE)
                    {
                        *l_pl_Wk++ = 0x0;
                        *l_pl_Jk_FG++ = 0x0;
                        *l_pl_Jk_BG2++ = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    }
                    else
                    {
                        *l_pl_Jk_BG2++ = 0;
                        *l_pl_Wk = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                        l_pl_Wk++;
                    }
#else
                    *l_pl_Wk = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                    l_pl_Wk++;
#endif
                    l_pt_D_DS++;
                    l_pt_Q_DS++;
                }
            }
        }
        {
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#ifdef	MLF_BackGround_Counter
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#endif
        }
        {
            const s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;

            s32* l_pus_J = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pus_J_BG = l_pl_Jk_BG;
#endif
            y = (s32)m_DataFusion.m_h*(s32)m_DataFusion.m_w;
            while (y--)
            {
                {
#ifdef	MLF_BackGround_Counter
                    if (*l_pl_Wk <= *l_pus_J_BG) // FGa, in the case where the foreground weight is 0, we set the background mask to 1
                        *l_pus_J_BG = 255;
                    else
                        *l_pus_J_BG = (255*(*l_pus_J_BG))/(2*(*l_pl_Wk));
                    l_pus_J_BG++;
#endif
                    if (*l_pl_Wk)
                    {
                        s32 l_lWkHalf = *l_pl_Wk/2; // FGa, to avoid numeric problems (in the test was between 999 and 1000)
                        *l_pus_J = ((*l_pus_J+l_lWkHalf))/(*l_pl_Wk);
                    }
                    l_pus_J++;
                    l_pl_Wk++;
                }
            }
        }
        // INTERPOLATION: Check which pixels in the high-resolution 2D image have an intensity value that belongs to the range (levels 0-1) that we've processed
        //                in order to update the final pixels
        if (i>0)
        {
            const s32 l_lXMax = (s32)((m_DataFusion.m_w-1)*l_ucShiftUnit);
            const s32 l_lYMax = (s32)((m_DataFusion.m_h-1)*l_ucShiftUnit)-1;
#ifdef	MLF_BackGround_Counter
            const s32* p_plImgJk0_BG = l_pt_Jk_BG[jk_0];
            const s32* p_plImgJk0_d_BG = p_plImgJk0_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_BG = l_pt_Jk_BG[jk_1];
            const s32* p_plImgJk1_d_BG = p_plImgJk1_BG+m_DataFusion.m_w;
#endif
            const s32* p_plImgJk0 = l_pt_Jk[jk_0];
            const s32* p_plImgJk0_d = p_plImgJk0+m_DataFusion.m_w;
            const s32* p_plImgJk1 = l_pt_Jk[jk_1];
            const s32* p_plImgJk1_d = p_plImgJk1+m_DataFusion.m_w;
            const s16* l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[0][0][0];
            const s32 l_lImgRangeXSihiftUnitSq = (s32)l_tImg_Range* (s32)l_ucShiftUnit*(s32)l_ucShiftUnit;
            const u08 *l_pt_I = m_DataFusion.m_pt_I_GRAY;
            const u16 l_tKf1Level = (((i-1)*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            const u16 l_tKf2Level = (((i  )*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            s16 l_sY0 = 0;
            s16 l_sdY = 0;
            u16* l_pt_J = m_DataFusion.m_pt_J;

            for (y=0;y<m_DataFusion.m_h_original;y++)
            {
                {
                    s16 l_sdX = 0;
                    u08 l_ucVal;

                    for (x=0;x<l_lXMax;x++)
                    {
                        l_ucVal = *l_pt_I;
                        if (l_ucVal<l_tKf2Level)
                        {
                           if (l_ucVal>=l_tKf1Level)
                           {
                                // qx_linear_interpolate_xy
                                s32 l_lRetInterpolate2;
                                s32 l_lRetInterpolate;
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d_BG+1));
                                l_lRetInterpolate *= alpha;
                                l_psaDxx -= 3; // reset the position
                                l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate2 *= alpha_Comp;
                                l_psaDxx -= 3; // reset the position
                                l_BakGroundWeight =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d+1));
                                    l_lRetInterpolate *= alpha;
                                    l_psaDxx -= 3;// reset the position
                                    l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate2 *= alpha_Comp;

                                    *l_pt_J =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                }
                            }
                        }
                        else
                        {
                            if (i==l_ucQuantizationLevels)
                            {
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
                                s32 l_lRetInterpolate ;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                l_psaDxx -= 3;// reset the position
                                l_BakGroundWeight = l_lRetInterpolate;
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                    l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                    *l_pt_J =(u16)(l_lRetInterpolate);
                                }
                            }
                        }
                        l_pt_J++;
                        l_pt_I++;
                        {
                            l_sdX++;
                            if (l_sdX == l_ucShiftUnit)
                            {
                                l_sdX= 0;
#ifdef	MLF_BackGround_Counter
                                p_plImgJk0_BG++;
                                p_plImgJk0_d_BG++;
                                p_plImgJk1_BG++;
                                p_plImgJk1_d_BG++;
#endif

                                p_plImgJk0++;
                                p_plImgJk0_d++;
                                p_plImgJk1++;
                                p_plImgJk1_d++;
                            }
                        }
                    }
                    // X Max to end
                    {
                        s32 l_dLIc;
                        s32 l_dLIa;
                        s32 l_dLIb;
                        const s16* l_psaDxx = l_spaDxx_Ylevel;
                        /*qx_linear_interpolate_xy */
#ifdef	MLF_BackGround_Counter
                        s32 l_dLIc_BG;
                        s32 l_dLIa_BG;
                        s32 l_dLIb_BG;
                        l_dLIa_BG = ((*l_psaDxx)  *(*(p_plImgJk0_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d_BG)));
                        l_dLIb_BG = (  (*l_psaDxx)*(*(p_plImgJk1_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d_BG)));
                        l_dLIc_BG = l_dLIb_BG+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
#endif
                        l_dLIa = ((*l_psaDxx)  *(*(p_plImgJk0))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d)));
                        l_dLIb = (  (*l_psaDxx)*(*(p_plImgJk1))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d)));
                        l_dLIc = l_dLIb+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
                        for (x = l_lXMax;x<m_DataFusion.m_w_original;x++)
                        {
                            if (l_ucVal<l_tKf2Level)
                            {
                              if (l_ucVal>=l_tKf1Level)
                              {
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);

#ifdef	MLF_BackGround_Counter
                                *l_pt_J =(u16)(((alpha*l_dLIa_BG+(alpha_Comp)* l_dLIb_BG))/(l_lImgRangeXSihiftUnitSq));
                                if (*l_pt_J > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                    *l_pt_J =(u16)(((alpha*l_dLIa+(alpha_Comp)* l_dLIb))/(l_lImgRangeXSihiftUnitSq));
                              }
                            }
                            else
                            {
                                if (i==l_ucQuantizationLevels)
                                {
#ifdef	MLF_BackGround_Counter
                                    *l_pt_J =(u16)(l_dLIc_BG/(l_ucShiftUnit*l_ucShiftUnit));
                                    if (*l_pt_J > D_BG_THRESHOLD)
                                        *l_pt_J = MLF_APP_MAX_DISTANCE;
                                    else
#endif
                                        *l_pt_J =(u16)(l_dLIc/(l_ucShiftUnit*l_ucShiftUnit));
                                }
                            }
                            l_pt_J++;
                            l_pt_I++;
                        }
                    }
                }
                if (y < l_lYMax)
                {
                    l_sdY++;
                    if (l_sdY == l_ucShiftUnit)
                    {
                        l_sdY= 0;
                        l_sY0++;
                    }
                    l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[l_sdY][0][0];
                }
#ifdef	MLF_BackGround_Counter
                p_plImgJk0_BG   = l_pt_Jk_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_BG = p_plImgJk0_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_BG   = l_pt_Jk_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_BG = p_plImgJk1_BG + (s32)m_DataFusion.m_w;
#endif

                p_plImgJk0   = l_pt_Jk[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d = p_plImgJk0 + (s32)m_DataFusion.m_w;
                p_plImgJk1   = l_pt_Jk[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d = p_plImgJk1 + (s32)m_DataFusion.m_w;
            }
            jk_1=jk_0;
            jk_0=(jk_0+1)%2;
        }
    }
#ifdef MLF_ANALYSE_TIME_CONSUMPTION
    }
    std::clock_t c_end = std::clock();

    std::cout << "CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC
              << " ms" << std::endl;
    std::cout << "Done..." << std::endl;
#endif

    return 0;
}

/* Apply Pixel Weighted Average Strategy RGB Filter (PWAS RGB Filter)
 */
short c_DataFusion::PWAS_RGB_Filter()
{
    const u16* l_pt_D_DS;
    const u08* l_pt_I_R_DS;
    const u08* l_pt_I_G_DS;
    const u08* l_pt_I_B_DS;
    const u08* l_pt_Q_DS;
    u08 l_tImage_min;
    u08 l_tImage_max;
    u08 l_tImg_Range;
    const u08 l_ucQuantizationLevels = D_QUANTIZATION_RANGE_2D-1;

    if (!m_DataFusion.m_sInitOk)
    {
        return -1;
    }
    if ((u16)D_QUANTIZATION_RANGE_2D > MLF_APP_INTENSITY_RANGE)
    {
        return -2;
    }
#ifdef	MLF_ANALYSE_TIME_CONSUMPTION
    std::cout << "Time consumption analysis for the PWAS RGB filter: " << std::endl;
    const clock_t c_start = std::clock();
    short l_sNbIter = 1;
    for (int ii=0;ii<l_sNbIter;ii++)
    {
#endif        
    // Compute the range [min,max] of the 2D image. This range will be divided into D_QUANTIZATION_RANGE_2D levels
    GetMinMaxVal_RGB(&l_tImage_min, &l_tImage_max, m_DataFusion.m_pt_I_R_DS, m_DataFusion.m_pt_I_G_DS, m_DataFusion.m_pt_I_B_DS, m_DataFusion.m_h*m_DataFusion.m_w);
    l_tImg_Range = l_tImage_max-l_tImage_min;
    s32 *l_pt_Jk_R[2];
    s32 *l_pl_Jk_R;
    s32 *l_pt_Jk_G[2];
    s32 *l_pl_Jk_G;
    s32 *l_pt_Jk_B[2];
    s32 *l_pl_Jk_B;
    s16 jk_0 = 0;
    s16 jk_1 = 1;
    s32 y,x;
    u08  l_ucShiftUnit = (1<<m_DataFusion.m_nr_shift);  // 1..16
#ifdef	MLF_BackGround_Counter
    s32 *l_pt_Jk_R_BG[2];
    s32 *l_pl_Jk_R_BG;
    l_pt_Jk_R_BG[0] = m_DataFusion.m_pt_Jk_level0_R_BG;
    l_pt_Jk_R_BG[1] = m_DataFusion.m_pt_Jk_level1_R_BG;
    s32 *l_pt_Jk_G_BG[2];
    s32 *l_pl_Jk_G_BG;
    l_pt_Jk_G_BG[0] = m_DataFusion.m_pt_Jk_level0_G_BG;
    l_pt_Jk_G_BG[1] = m_DataFusion.m_pt_Jk_level1_G_BG;
    s32 *l_pt_Jk_B_BG[2];
    s32 *l_pl_Jk_B_BG;
    l_pt_Jk_B_BG[0] = m_DataFusion.m_pt_Jk_level0_B_BG;
    l_pt_Jk_B_BG[1] = m_DataFusion.m_pt_Jk_level1_B_BG;
#endif
    l_pt_Jk_R[0] = m_DataFusion.m_pt_Jk_level0_R;
    l_pt_Jk_R[1] = m_DataFusion.m_pt_Jk_level1_R;
    l_pt_Jk_G[0] = m_DataFusion.m_pt_Jk_level0_G;
    l_pt_Jk_G[1] = m_DataFusion.m_pt_Jk_level1_G;
    l_pt_Jk_B[0] = m_DataFusion.m_pt_Jk_level0_B;
    l_pt_Jk_B[1] = m_DataFusion.m_pt_Jk_level1_B;

    // For all levels
    for (s16 i=0;i<D_QUANTIZATION_RANGE_2D;i++)
    {
        if (i==0)
        {
            l_pl_Jk_R = l_pt_Jk_R[jk_0];
            l_pl_Jk_G = l_pt_Jk_G[jk_0];
            l_pl_Jk_B = l_pt_Jk_B[jk_0];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_R_BG = l_pt_Jk_R_BG[jk_0];
            l_pl_Jk_G_BG = l_pt_Jk_G_BG[jk_0];
            l_pl_Jk_B_BG = l_pt_Jk_B_BG[jk_0];
#endif
        }
        else
        {
            l_pl_Jk_R = l_pt_Jk_R[jk_1];
            l_pl_Jk_G = l_pt_Jk_G[jk_1];
            l_pl_Jk_B = l_pt_Jk_B[jk_1];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_R_BG = l_pt_Jk_R_BG[jk_1];
            l_pl_Jk_G_BG = l_pt_Jk_G_BG[jk_1];
            l_pl_Jk_B_BG = l_pt_Jk_B_BG[jk_1];
#endif
        }
        {
            s32* l_pl_Wk_R = m_DataFusion.m_pt_Wk_R;
            s32* l_pl_Wk_G = m_DataFusion.m_pt_Wk_G;
            s32* l_pl_Wk_B = m_DataFusion.m_pt_Wk_B;
            s32* l_pl_Jk_R_FG = l_pl_Jk_R;
            s32* l_pl_Jk_G_FG = l_pl_Jk_G;
            s32* l_pl_Jk_B_FG = l_pl_Jk_B;
#ifdef	MLF_BackGround_Counter
            s32* l_pl_Jk_R_BG2 = l_pl_Jk_R_BG;
            s32* l_pl_Jk_G_BG2 = l_pl_Jk_G_BG;
            s32* l_pl_Jk_B_BG2 = l_pl_Jk_B_BG;
#endif
            // RGB to Grayscale
            s32 l_lIndexIni;
            if (i== 0)
                l_lIndexIni = l_tImage_min;
            else
            {
                if (i == l_ucQuantizationLevels)
                    l_lIndexIni = l_tImage_max;
                else
                    l_lIndexIni = ((s32)l_tImage_min +(i*(s32)(l_tImg_Range))/l_ucQuantizationLevels);
            }
            l_pt_D_DS = m_DataFusion.m_pt_D_DS;
            l_pt_I_R_DS = m_DataFusion.m_pt_I_R_DS;
            l_pt_I_G_DS = m_DataFusion.m_pt_I_G_DS;
            l_pt_I_B_DS = m_DataFusion.m_pt_I_B_DS;
            l_pt_Q_DS = m_DataFusion.m_pt_Q_DS;

            for (y=0;y<m_DataFusion.m_h;y++)
            {
                for (x=0;x<m_DataFusion.m_w;x++)
                {
                    s32 index_R=(l_lIndexIni-(s32)(*l_pt_I_R_DS++));
                    s32 index_G=(l_lIndexIni-(s32)(*l_pt_I_G_DS++));
                    s32 index_B=(l_lIndexIni-(s32)(*l_pt_I_B_DS++));
                    if(index_R < 0)
                        index_R = -index_R;
                    if(index_G < 0)
                        index_G = -index_G;
                    if(index_B < 0)
                        index_B = -index_B;
#ifdef	MLF_BackGround_Counter
                    if (*l_pt_D_DS >= MLF_APP_MAX_DISTANCE)
                    {
                        *l_pl_Wk_R++ = 0x0;
                        *l_pl_Jk_R_FG++ = 0x0;
                        *l_pl_Jk_R_BG2++ = (s32)(m_DataFusion.m_ucaTable[index_R]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Wk_G++ = 0x0;
                        *l_pl_Jk_G_FG++ = 0x0;
                        *l_pl_Jk_G_BG2++ = (s32)(m_DataFusion.m_ucaTable[index_G]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Wk_B++ = 0x0;
                        *l_pl_Jk_B_FG++ = 0x0;
                        *l_pl_Jk_B_BG2++ = (s32)(m_DataFusion.m_ucaTable[index_B]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    }
                    else
                    {
                        *l_pl_Jk_R_BG2++ = 0x0;
                        *l_pl_Wk_R = (s32)(m_DataFusion.m_ucaTable[index_R]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_R_FG++ = (*l_pl_Wk_R)*(*l_pt_D_DS);
                        l_pl_Wk_R++;
                        *l_pl_Jk_G_BG2++ = 0x0;
                        *l_pl_Wk_G = (s32)(m_DataFusion.m_ucaTable[index_G]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_G_FG++ = (*l_pl_Wk_G)*(*l_pt_D_DS);
                        l_pl_Wk_G++;
                        *l_pl_Jk_B_BG2++ = 0x0;
                        *l_pl_Wk_B = (s32)(m_DataFusion.m_ucaTable[index_B]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_B_FG++ = (*l_pl_Wk_B)*(*l_pt_D_DS);
                        l_pl_Wk_B++;
                    }
#else
                    *l_pl_Wk_R = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_R_FG++ = (*l_pl_Wk_R)*(*l_pt_D_DS);
                    l_pl_Wk_R++;
                    *l_pl_Wk_G = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_G_FG++ = (*l_pl_Wk_G)*(*l_pt_D_DS);
                    l_pl_Wk_G++;
                    *l_pl_Wk_B = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_B_FG++ = (*l_pl_Wk_B)*(*l_pt_D_DS);
                    l_pl_Wk_B++;
#endif
                    l_pt_D_DS++;
                    l_pt_Q_DS++;
                }
            }
        }
        {
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_G,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_G,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_B,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_B,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#ifdef	MLF_BackGround_Counter
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_R_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_G_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_B_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#endif
        }
        {
            const s32* l_pl_Wk_R = m_DataFusion.m_pt_Wk_R;
            const s32* l_pl_Wk_G = m_DataFusion.m_pt_Wk_G;
            const s32* l_pl_Wk_B = m_DataFusion.m_pt_Wk_B;
            s32* l_pus_J_R = l_pl_Jk_R;
            s32* l_pus_J_G = l_pl_Jk_G;
            s32* l_pus_J_B = l_pl_Jk_B;
#ifdef	MLF_BackGround_Counter
            s32* l_pus_J_R_BG = l_pl_Jk_R_BG;
            s32* l_pus_J_G_BG = l_pl_Jk_G_BG;
            s32* l_pus_J_B_BG = l_pl_Jk_B_BG;
#endif
            y = (s32)m_DataFusion.m_h*(s32)m_DataFusion.m_w;
            while (y--)
            {
                {
#ifdef	MLF_BackGround_Counter
                    if (*l_pl_Wk_R <= *l_pus_J_R_BG) // FGa, in the case where the foreground weight is 0, we set the background mask to 1
                        *l_pus_J_R_BG = 255;
                    else
                        *l_pus_J_R_BG = (255*(*l_pus_J_R_BG))/(2*(*l_pl_Wk_R));
                    if (*l_pl_Wk_G <= *l_pus_J_G_BG)
                        *l_pus_J_G_BG = 255;
                    else
                        *l_pus_J_G_BG = (255*(*l_pus_J_G_BG))/(2*(*l_pl_Wk_G));
                    if (*l_pl_Wk_B <= *l_pus_J_B_BG)
                        *l_pus_J_B_BG = 255;
                    else
                        *l_pus_J_B_BG = (255*(*l_pus_J_B_BG))/(2*(*l_pl_Wk_B));
                    l_pus_J_R_BG++;
                    l_pus_J_G_BG++;
                    l_pus_J_B_BG++;
#endif
                    if (*l_pl_Wk_R)
                        *l_pus_J_R = ((*l_pus_J_R))/(*l_pl_Wk_R);
                    if (*l_pl_Wk_G)
                        *l_pus_J_G = ((*l_pus_J_G))/(*l_pl_Wk_G);
                    if (*l_pl_Wk_B)
                        *l_pus_J_B = ((*l_pus_J_B))/(*l_pl_Wk_B);
                    l_pl_Wk_R++;
                    l_pl_Wk_G++;
                    l_pl_Wk_B++;
                    l_pus_J_R++;
                    l_pus_J_G++;
                    l_pus_J_B++;
                }
            }
        }
        // INTERPOLATION: Check which pixels in the high-resolution 2D image have an intensity value that belongs to the range (levels 0-1) that we've processed
        //                in order to update the final pixels
        if (i>0)
        {
            const s32 l_lXMax = (s32)((m_DataFusion.m_w-1)*l_ucShiftUnit);
            const s32 l_lYMax = (s32)((m_DataFusion.m_h-1)*l_ucShiftUnit)-1;
#ifdef	MLF_BackGround_Counter
            const s32* p_plImgJk0_BG;
            const s32* p_plImgJk0_R_BG = l_pt_Jk_R_BG[jk_0];
            const s32* p_plImgJk0_G_BG = l_pt_Jk_G_BG[jk_0];
            const s32* p_plImgJk0_B_BG = l_pt_Jk_B_BG[jk_0];
            const s32* p_plImgJk0_d_BG;
            const s32* p_plImgJk0_d_R_BG = p_plImgJk0_R_BG+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_G_BG = p_plImgJk0_G_BG+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_B_BG = p_plImgJk0_B_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_BG;
            const s32* p_plImgJk1_R_BG = l_pt_Jk_R_BG[jk_1];
            const s32* p_plImgJk1_G_BG = l_pt_Jk_G_BG[jk_1];
            const s32* p_plImgJk1_B_BG = l_pt_Jk_B_BG[jk_1];
            const s32* p_plImgJk1_d_BG;
            const s32* p_plImgJk1_d_R_BG = p_plImgJk1_R_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_G_BG = p_plImgJk1_G_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_B_BG = p_plImgJk1_B_BG+m_DataFusion.m_w;
#endif
            const s32* p_plImgJk0;
            const s32* p_plImgJk0_R = l_pt_Jk_R[jk_0];
            const s32* p_plImgJk0_G = l_pt_Jk_G[jk_0];
            const s32* p_plImgJk0_B = l_pt_Jk_B[jk_0];
            const s32* p_plImgJk0_d;
            const s32* p_plImgJk0_d_R = p_plImgJk0_R+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_G = p_plImgJk0_G+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_B = p_plImgJk0_B+m_DataFusion.m_w;
            const s32* p_plImgJk1;
            const s32* p_plImgJk1_R = l_pt_Jk_R[jk_1];
            const s32* p_plImgJk1_G = l_pt_Jk_G[jk_1];
            const s32* p_plImgJk1_B = l_pt_Jk_B[jk_1];
            const s32* p_plImgJk1_d;
            const s32* p_plImgJk1_d_R = p_plImgJk1_R+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_G = p_plImgJk1_G+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_B = p_plImgJk1_B+m_DataFusion.m_w;
            const s16* l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[0][0][0];
            const s32 l_lImgRangeXSihiftUnitSq = (s32)l_tImg_Range* (s32)l_ucShiftUnit*(s32)l_ucShiftUnit;
            const u08 *l_pt_I_R = m_DataFusion.m_pt_I_R;
            const u08 *l_pt_I_G = m_DataFusion.m_pt_I_G;
            const u08 *l_pt_I_B = m_DataFusion.m_pt_I_B;
            const u08 *l_pt_BetaChannel = m_DataFusion.m_pt_BetaChannel;
            const u16 l_tKf1Level = (((i-1)*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            const u16 l_tKf2Level = (((i  )*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            s16 l_sY0 = 0;
            s16 l_sdY = 0;
            u16* l_pt_J = m_DataFusion.m_pt_J;

            for (y=0;y<m_DataFusion.m_h_original;y++)
            {
                {
                    s16 l_sdX = 0;
                    u08 l_ucVal;

                    for (x=0;x<l_lXMax;x++)
                    {
                        if (*l_pt_BetaChannel == D_CHANNEL_R) // R-Channel
                        {
                            l_ucVal = *l_pt_I_R;
                            p_plImgJk0_BG = p_plImgJk0_R_BG;
                            p_plImgJk0_d_BG = p_plImgJk0_d_R_BG;
                            p_plImgJk1_BG = p_plImgJk1_R_BG;
                            p_plImgJk1_d_BG = p_plImgJk1_d_R_BG;
                            p_plImgJk0 = p_plImgJk0_R;
                            p_plImgJk0_d = p_plImgJk0_d_R;
                            p_plImgJk1 = p_plImgJk1_R;
                            p_plImgJk1_d = p_plImgJk1_d_R;

                        }
                        else if (*l_pt_BetaChannel == D_CHANNEL_G) // G-Channel
                        {
                            l_ucVal = *l_pt_I_G;
                            p_plImgJk0_BG = p_plImgJk0_G_BG;
                            p_plImgJk0_d_BG = p_plImgJk0_d_G_BG;
                            p_plImgJk1_BG = p_plImgJk1_G_BG;
                            p_plImgJk1_d_BG = p_plImgJk1_d_G_BG;
                            p_plImgJk0 = p_plImgJk0_G;
                            p_plImgJk0_d = p_plImgJk0_d_G;
                            p_plImgJk1 = p_plImgJk1_G;
                            p_plImgJk1_d = p_plImgJk1_d_G;
                        }
                        else if (*l_pt_BetaChannel == D_CHANNEL_B) // B-Channel
                        {
                            l_ucVal = *l_pt_I_B;
                            p_plImgJk0_BG = p_plImgJk0_B_BG;
                            p_plImgJk0_d_BG = p_plImgJk0_d_B_BG;
                            p_plImgJk1_BG = p_plImgJk1_B_BG;
                            p_plImgJk1_d_BG = p_plImgJk1_d_B_BG;
                            p_plImgJk0 = p_plImgJk0_B;
                            p_plImgJk0_d = p_plImgJk0_d_B;
                            p_plImgJk1 = p_plImgJk1_B;
                            p_plImgJk1_d = p_plImgJk1_d_B;
                        }
                        if (l_ucVal<l_tKf2Level)
                        {
                           if (l_ucVal>=l_tKf1Level)
                           {
                                // qx_linear_interpolate_xy
                                s32 l_lRetInterpolate2;
                                s32 l_lRetInterpolate;
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d_BG+1));
                                l_lRetInterpolate *= alpha;
                                l_psaDxx -= 3; // reset the position
                                l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate2 *= alpha_Comp;
                                l_psaDxx -= 3; // reset the position
                                l_BakGroundWeight =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d+1));
                                    l_lRetInterpolate *= alpha;

                                    l_psaDxx -= 3;// reset the position
                                    l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate2 *= alpha_Comp;

                                    *l_pt_J =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                }
                            }
                        }
                        else
                        {
                            if (i==l_ucQuantizationLevels)
                            {
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
                                s32 l_lRetInterpolate ;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                l_psaDxx -= 3;// reset the position
                                l_BakGroundWeight = l_lRetInterpolate;
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                    l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                    *l_pt_J =(u16)(l_lRetInterpolate);
                                }
                            }
                        }
                        l_pt_J++;
                        l_pt_I_R++;
                        l_pt_I_B++;
                        l_pt_I_G++;
                        l_pt_BetaChannel++;
                        {
                            l_sdX++;
                            if (l_sdX == l_ucShiftUnit)
                            {
                                l_sdX= 0;
#ifdef	MLF_BackGround_Counter
                                p_plImgJk0_R_BG++;
                                p_plImgJk0_G_BG++;
                                p_plImgJk0_B_BG++;
                                p_plImgJk0_d_R_BG++;
                                p_plImgJk0_d_G_BG++;
                                p_plImgJk0_d_B_BG++;
                                p_plImgJk1_R_BG++;
                                p_plImgJk1_G_BG++;
                                p_plImgJk1_B_BG++;
                                p_plImgJk1_d_R_BG++;
                                p_plImgJk1_d_G_BG++;
                                p_plImgJk1_d_B_BG++;
#endif

                                p_plImgJk0_R++;
                                p_plImgJk0_G++;
                                p_plImgJk0_B++;
                                p_plImgJk0_d_R++;
                                p_plImgJk0_d_G++;
                                p_plImgJk0_d_B++;
                                p_plImgJk1_R++;
                                p_plImgJk1_G++;
                                p_plImgJk1_B++;
                                p_plImgJk1_d_R++;
                                p_plImgJk1_d_G++;
                                p_plImgJk1_d_B++;
                            }
                        }
                    }
                    // X Max to end
                    {
                        s32 l_dLIc;
                        s32 l_dLIa;
                        s32 l_dLIb;
                        const s16* l_psaDxx = l_spaDxx_Ylevel;
                        /*qx_linear_interpolate_xy */
#ifdef	MLF_BackGround_Counter
                        s32 l_dLIc_BG;
                        s32 l_dLIa_BG;
                        s32 l_dLIb_BG;
                        l_dLIa_BG = ((*l_psaDxx)  *(*(p_plImgJk0_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d_BG)));
                        l_dLIb_BG = (  (*l_psaDxx)*(*(p_plImgJk1_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d_BG)));
                        l_dLIc_BG = l_dLIb_BG+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
#endif
                        l_dLIa = ((*l_psaDxx)  *(*(p_plImgJk0))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d)));
                        l_dLIb = (  (*l_psaDxx)*(*(p_plImgJk1))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d)));
                        l_dLIc = l_dLIb+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
                        for (x = l_lXMax;x<m_DataFusion.m_w_original;x++)
                        {
                            if (l_ucVal<l_tKf2Level)
                            {
                              if (l_ucVal>=l_tKf1Level)
                              {
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);

#ifdef	MLF_BackGround_Counter
                                *l_pt_J =(u16)(((alpha*l_dLIa_BG+(alpha_Comp)* l_dLIb_BG))/(l_lImgRangeXSihiftUnitSq));
                                if (*l_pt_J > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                    *l_pt_J =(u16)(((alpha*l_dLIa+(alpha_Comp)* l_dLIb))/(l_lImgRangeXSihiftUnitSq));
                              }
                            }
                            else
                            {
                                if (i==l_ucQuantizationLevels)
                                {
#ifdef	MLF_BackGround_Counter
                                    *l_pt_J =(u16)(l_dLIc_BG/(l_ucShiftUnit*l_ucShiftUnit));
                                    if (*l_pt_J > D_BG_THRESHOLD)
                                        *l_pt_J = MLF_APP_MAX_DISTANCE;
                                    else
#endif
                                        *l_pt_J =(u16)(l_dLIc/(l_ucShiftUnit*l_ucShiftUnit));
                                }
                            }
                            l_pt_J++;
                            l_pt_I_R++;
                            l_pt_I_B++;
                            l_pt_I_G++;
                            l_pt_BetaChannel++;
                        }
                    }
                }
                if (y < l_lYMax)
                {
                    l_sdY++;
                    if (l_sdY == l_ucShiftUnit)
                    {
                        l_sdY= 0;
                        l_sY0++;
                    }
                    l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[l_sdY][0][0];
                }
#ifdef	MLF_BackGround_Counter
                p_plImgJk0_R_BG   = l_pt_Jk_R_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_G_BG   = l_pt_Jk_G_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_B_BG   = l_pt_Jk_B_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_R_BG = p_plImgJk0_R_BG + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_G_BG = p_plImgJk0_G_BG + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_B_BG = p_plImgJk0_B_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_R_BG   = l_pt_Jk_R_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_G_BG   = l_pt_Jk_G_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_B_BG   = l_pt_Jk_B_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_R_BG = p_plImgJk1_R_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_G_BG = p_plImgJk1_G_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_B_BG = p_plImgJk1_B_BG + (s32)m_DataFusion.m_w;
#endif

                p_plImgJk0_R   = l_pt_Jk_R[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_G   = l_pt_Jk_G[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_B   = l_pt_Jk_B[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_R = p_plImgJk0_R + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_G = p_plImgJk0_G + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_B = p_plImgJk0_B + (s32)m_DataFusion.m_w;
                p_plImgJk1_R   = l_pt_Jk_R[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_G   = l_pt_Jk_G[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_B   = l_pt_Jk_B[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_R = p_plImgJk1_R + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_G = p_plImgJk1_G + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_B = p_plImgJk1_B + (s32)m_DataFusion.m_w;
            }
            jk_1=jk_0;
            jk_0=(jk_0+1)%2;
        }
    }
#ifdef MLF_ANALYSE_TIME_CONSUMPTION
    }
    std::clock_t c_end = std::clock();

    std::cout << "CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC
              << " ms" << std::endl;
    std::cout << "Done..." << std::endl;
#endif

    return 0;
}

/* Apply Unified Multi-Lateral Filter (UML Filter)
 */
short c_DataFusion::UML_Filter()
{
    const u16* l_pt_D_DS;
    const u08* l_pt_I_DS;
    const u08* l_pt_Q_DS;
    u08 l_tImage_min;
    u08 l_tImage_max;
    u08 l_tImg_Range;
    const u08 l_ucQuantizationLevels = D_QUANTIZATION_RANGE_2D-1;

    if (!m_DataFusion.m_sInitOk)
    {
        return -1;
    }
    if ((u16)D_QUANTIZATION_RANGE_2D > MLF_APP_2D_RANGE)
    {
        return -2;
    }
#ifdef	MLF_ANALYSE_TIME_CONSUMPTION
    std::cout << "Time consumption analysis for the UML filter: " << std::endl;
    const clock_t c_start = std::clock();
    short l_sNbIter = 1;
    for (int ii=0;ii<l_sNbIter;ii++)
    {
#endif
    // Compute the range [min,max] of the 2D image. This range will be divided into D_QUANTIZATION_RANGE_2D levels
    GetMinMaxVal(&l_tImage_min, &l_tImage_max, m_DataFusion.m_pt_I_GRAY_DS, m_DataFusion.m_h*m_DataFusion.m_w);
    l_tImg_Range = l_tImage_max-l_tImage_min;
    s32 *l_pt_Jk[2];
    s32 *l_pl_Jk;
    s16 jk_0 = 0;
    s16 jk_1 = 1;
    s32 y,x;
    u08  l_ucShiftUnit = (1<<m_DataFusion.m_nr_shift);  // 1..16
#ifdef	MLF_BackGround_Counter
    s32 *l_pt_Jk_BG[2];
    s32 *l_pl_Jk_BG;
    l_pt_Jk_BG[0] = m_DataFusion.m_pt_Jk_level0_R_BG;
    l_pt_Jk_BG[1] = m_DataFusion.m_pt_Jk_level1_R_BG;
#endif
    l_pt_Jk[0] = m_DataFusion.m_pt_Jk_level0_R;
    l_pt_Jk[1] = m_DataFusion.m_pt_Jk_level1_R;

    // For all levels
    for (s16 i=0;i<D_QUANTIZATION_RANGE_2D;i++)
    {
        if (i==0)
        {
            l_pl_Jk=	l_pt_Jk[jk_0];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_0];
#endif
        }
        else
        {
            l_pl_Jk=	l_pt_Jk[jk_1];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_BG=l_pt_Jk_BG[jk_1];
#endif
        }
        {
            s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;
            s32* l_pl_Jk_FG = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pl_Jk_BG2 = l_pl_Jk_BG;
#endif
            // RGB to Grayscale
            s32 l_lIndexIni;
            if (i== 0)
                l_lIndexIni = l_tImage_min;
            else
            {
                if (i == l_ucQuantizationLevels)
                    l_lIndexIni = l_tImage_max;
                else
                    l_lIndexIni = ((s32)l_tImage_min +(i*(s32)(l_tImg_Range))/l_ucQuantizationLevels);
            }
            l_pt_D_DS = m_DataFusion.m_pt_D_DS;
            l_pt_I_DS = m_DataFusion.m_pt_I_GRAY_DS;
            l_pt_Q_DS = m_DataFusion.m_pt_Q_DS;

            for (y=0;y<m_DataFusion.m_h;y++)
            {
                for (x=0;x<m_DataFusion.m_w;x++)
                {
                    s32 index=(l_lIndexIni-(s32)(*l_pt_I_DS++));
                    if (index < 0)
                    {
                        index = -index;
                    }/*
                    if (index >= MLF_APP_2D_RANGE) // Updated by BMi
                    {
                        index = MLF_APP_2D_RANGE-1;
                    }*/
#ifdef	MLF_BackGround_Counter
                    if (*l_pt_D_DS >= MLF_APP_MAX_DISTANCE)
                    {
                        *l_pl_Wk++ = 0x0;
                        *l_pl_Jk_FG++ = 0x0;
                        *l_pl_Jk_BG2++ = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    }
                    else
                    {
                        *l_pl_Jk_BG2++ = 0;
                        *l_pl_Wk = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                        l_pl_Wk++;
                    }
#else
                    *l_pl_Wk = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_FG++ = (*l_pl_Wk)*(*l_pt_D_DS);
                    l_pl_Wk++;
#endif
                    l_pt_D_DS++;
                    l_pt_Q_DS++;
                }
            }
        }
        {
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#ifdef	MLF_BackGround_Counter
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#endif
        }
        {
            const s32* l_pl_Wk = m_DataFusion.m_pt_Wk_R;

            s32* l_pus_J = l_pl_Jk;
#ifdef	MLF_BackGround_Counter
            s32* l_pus_J_BG = l_pl_Jk_BG;
#endif
            y = (s32)m_DataFusion.m_h*(s32)m_DataFusion.m_w;
            while (y--)
            {
                {
#ifdef	MLF_BackGround_Counter
                    if (*l_pl_Wk <= *l_pus_J_BG) // FGa, in the case where the foreground weight is 0, we set the background mask to 1
                        *l_pus_J_BG = 255;
                    else
                        *l_pus_J_BG = (255*(*l_pus_J_BG))/(2*(*l_pl_Wk));
                    l_pus_J_BG++;
#endif
                    if (*l_pl_Wk)
                    {
                        s32 l_lWkHalf = *l_pl_Wk/2; // FGa, to avoid numeric problems (in the test was between 999 and 1000)
                        *l_pus_J = ((*l_pus_J+l_lWkHalf))/(*l_pl_Wk);
                    }
                    l_pus_J++;
                    l_pl_Wk++;
                }
            }
        }
        // INTERPOLATION: Check which pixels in the high-resolution 2D image have an intensity value that belongs to the range (levels 0-1) that we've processed
        //                in order to update the final pixels
        if (i>0)
        {
            const s32 l_lXMax = (s32)((m_DataFusion.m_w-1)*l_ucShiftUnit);
            const s32 l_lYMax = (s32)((m_DataFusion.m_h-1)*l_ucShiftUnit)-1;
#ifdef	MLF_BackGround_Counter
            const s32* p_plImgJk0_BG = l_pt_Jk_BG[jk_0];
            const s32* p_plImgJk0_d_BG = p_plImgJk0_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_BG = l_pt_Jk_BG[jk_1];
            const s32* p_plImgJk1_d_BG = p_plImgJk1_BG+m_DataFusion.m_w;
#endif
            const s32* p_plImgJk0 = l_pt_Jk[jk_0];
            const s32* p_plImgJk0_d = p_plImgJk0+m_DataFusion.m_w;
            const s32* p_plImgJk1 = l_pt_Jk[jk_1];
            const s32* p_plImgJk1_d = p_plImgJk1+m_DataFusion.m_w;
            const s16* l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[0][0][0];
            const s32 l_lImgRangeXSihiftUnitSq = (s32)l_tImg_Range* (s32)l_ucShiftUnit*(s32)l_ucShiftUnit;
            const u08 *l_pt_I = m_DataFusion.m_pt_I_GRAY;
            const u16 *l_pt_D = m_DataFusion.m_pt_D;
            const u08 *l_pt_BetaValue = m_DataFusion.m_pt_Q;
            const u16 l_tKf1Level = (((i-1)*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            const u16 l_tKf2Level = (((i  )*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            s16 l_sY0 = 0;
            s16 l_sdY = 0;
            u16* l_pt_J = m_DataFusion.m_pt_J;

            for (y=0;y<m_DataFusion.m_h_original;y++)
            {
                {
                    s16 l_sdX = 0;
                    u08 l_ucVal;

                    for (x=0;x<l_lXMax;x++)
                    {
                        l_ucVal = *l_pt_I;
                        if (l_ucVal<l_tKf2Level)
                        {
                           if (l_ucVal>=l_tKf1Level)
                           {
                                // qx_linear_interpolate_xy
                                s32 l_lRetInterpolate2;
                                s32 l_lRetInterpolate;
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++) *(*(p_plImgJk0_BG));
                                l_lRetInterpolate += (*l_psaDxx++) *(*(p_plImgJk0_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++) *(*(p_plImgJk0_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)   *(*(p_plImgJk0_d_BG+1));
                                l_lRetInterpolate *= alpha;
                                l_psaDxx -= 3; // reset the position
                                l_lRetInterpolate2  = (*l_psaDxx++) *(*(p_plImgJk1_BG));
                                l_lRetInterpolate2 += (*l_psaDxx++) *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate2 += (*l_psaDxx++) *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate2 += (*l_psaDxx)   *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate2 *= alpha_Comp;
                                l_psaDxx -= 3; // reset the position
                                l_BakGroundWeight =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));

                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                {
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                }
                                else
#endif
                                {
                                    l_lRetInterpolate = (*l_psaDxx++)  *(*(p_plImgJk0));
                                    l_lRetInterpolate += (*l_psaDxx++) *(*(p_plImgJk0+1));
                                    l_lRetInterpolate += (*l_psaDxx++) *(*(p_plImgJk0_d));
                                    l_lRetInterpolate += (*l_psaDxx)   *(*(p_plImgJk0_d+1));
                                    l_lRetInterpolate *= alpha;
                                    l_psaDxx -= 3;// reset the position
                                    l_lRetInterpolate2 = (*l_psaDxx++)  *(*(p_plImgJk1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)  *(*(p_plImgJk1+1));
                                    l_lRetInterpolate2 += (*l_psaDxx++) *(*(p_plImgJk1_d));
                                    l_lRetInterpolate2 += (*l_psaDxx) *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate2 *= alpha_Comp;

                                    if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                        *l_pt_J =((u16)(255-(*l_pt_BetaValue))*((l_lRetInterpolate+l_lRetInterpolate2)/(l_lImgRangeXSihiftUnitSq)) + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                    else
                                        *l_pt_J =(u16)((l_lRetInterpolate+l_lRetInterpolate2)/(l_lImgRangeXSihiftUnitSq));
                                }
                            }
                        }
                        else
                        {
                            if (i==l_ucQuantizationLevels)
                            {
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
                                s32 l_lRetInterpolate ;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                l_psaDxx -= 3;// reset the position
                                l_BakGroundWeight = l_lRetInterpolate;
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                    l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                    if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                        *l_pt_J =((u16)(255-(*l_pt_BetaValue))*l_lRetInterpolate + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                    else
                                        *l_pt_J =(u16)(l_lRetInterpolate);
                                }
                            }
                        }
                        l_pt_J++;
                        l_pt_I++;
                        l_pt_D++;
                        l_pt_BetaValue++;
                        {
                            l_sdX++;
                            if (l_sdX == l_ucShiftUnit)
                            {
                                l_sdX= 0;
#ifdef	MLF_BackGround_Counter
                                p_plImgJk0_BG++;
                                p_plImgJk0_d_BG++;
                                p_plImgJk1_BG++;
                                p_plImgJk1_d_BG++;
#endif

                                p_plImgJk0++;
                                p_plImgJk0_d++;
                                p_plImgJk1++;
                                p_plImgJk1_d++;
                            }
                        }
                    }
                    // X Max to end
                    {
                        s32 l_dLIc;
                        s32 l_dLIa;
                        s32 l_dLIb;
                        const s16* l_psaDxx = l_spaDxx_Ylevel;
                        /*qx_linear_interpolate_xy */
#ifdef	MLF_BackGround_Counter
                        s32 l_dLIc_BG;
                        s32 l_dLIa_BG;
                        s32 l_dLIb_BG;
                        l_dLIa_BG = ((*l_psaDxx)  *(*(p_plImgJk0_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d_BG)));
                        l_dLIb_BG = (  (*l_psaDxx)*(*(p_plImgJk1_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d_BG)));
                        l_dLIc_BG = l_dLIb_BG+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
#endif
                        l_dLIa = ((*l_psaDxx)  *(*(p_plImgJk0))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d)));
                        l_dLIb = (  (*l_psaDxx)*(*(p_plImgJk1))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d)));
                        l_dLIc = l_dLIb+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
                        for (x = l_lXMax;x<m_DataFusion.m_w_original;x++)
                        {
                            if (l_ucVal<l_tKf2Level)
                            {
                              if (l_ucVal>=l_tKf1Level)
                              {
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);

#ifdef	MLF_BackGround_Counter
                                *l_pt_J =(u16)(((alpha*l_dLIa_BG+(alpha_Comp)* l_dLIb_BG))/(l_lImgRangeXSihiftUnitSq));
                                if (*l_pt_J > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                        *l_pt_J =((u16)(255-(*l_pt_BetaValue)) * (((alpha*l_dLIa+(alpha_Comp)*l_dLIb))/(l_lImgRangeXSihiftUnitSq)) + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                    else
                                        *l_pt_J =(u16)(((alpha*l_dLIa+(alpha_Comp) * l_dLIb))/(l_lImgRangeXSihiftUnitSq));
                                }
                              }
                            }
                            else
                            {
                                if (i==l_ucQuantizationLevels)
                                {
#ifdef	MLF_BackGround_Counter
                                    *l_pt_J =(u16)(l_dLIc_BG/(l_ucShiftUnit*l_ucShiftUnit));
                                    if (*l_pt_J > D_BG_THRESHOLD)
                                        *l_pt_J = MLF_APP_MAX_DISTANCE;
                                    else
#endif
                                        if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                            *l_pt_J =((u16)(255-(*l_pt_BetaValue))* (l_dLIc/(l_ucShiftUnit*l_ucShiftUnit)) + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                        else
                                            *l_pt_J =(u16)(l_dLIc/(l_ucShiftUnit*l_ucShiftUnit));
                                }
                            }
                            l_pt_J++;
                            l_pt_I++;
                            l_pt_D++;
                            l_pt_BetaValue++;
                        }
                    }
                }
                if (y < l_lYMax)
                {
                    l_sdY++;
                    if (l_sdY == l_ucShiftUnit)
                    {
                        l_sdY= 0;
                        l_sY0++;
                    }
                    l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[l_sdY][0][0];
                }
#ifdef	MLF_BackGround_Counter
                p_plImgJk0_BG   = l_pt_Jk_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_BG = p_plImgJk0_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_BG   = l_pt_Jk_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_BG = p_plImgJk1_BG + (s32)m_DataFusion.m_w;
#endif

                p_plImgJk0   = l_pt_Jk[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d = p_plImgJk0 + (s32)m_DataFusion.m_w;
                p_plImgJk1   = l_pt_Jk[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d = p_plImgJk1 + (s32)m_DataFusion.m_w;
            }
            jk_1=jk_0;
            jk_0=(jk_0+1)%2;
        }
    }
#ifdef MLF_ANALYSE_TIME_CONSUMPTION
    }
    std::clock_t c_end = std::clock();

    std::cout << "CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC
              << " ms" << std::endl;
    std::cout << "Done..." << std::endl;
#endif

    return 0;
}

/* Apply Unified Multi-Lateral RGB Filter (UML RGB Filter)
 */
short c_DataFusion::UML_RGB_Filter()
{
    const u16* l_pt_D_DS;
    const u08* l_pt_I_R_DS;
    const u08* l_pt_I_G_DS;
    const u08* l_pt_I_B_DS;
    const u08* l_pt_I_GRAY_DS;
    const u08* l_pt_Q_DS;
    u08 l_tImage_min;
    u08 l_tImage_max;
    u08 l_tImg_Range;
    const u08 l_ucQuantizationLevels = D_QUANTIZATION_RANGE_2D-1;

    if (!m_DataFusion.m_sInitOk)
    {
        return -1;
    }
    if ((u16)D_QUANTIZATION_RANGE_2D > MLF_APP_INTENSITY_RANGE)
    {
        return -2;
    }
#ifdef	MLF_ANALYSE_TIME_CONSUMPTION
    std::cout << "Time consumption analysis for the UML RGB filter: " << std::endl;
    const clock_t c_start = std::clock();
    short l_sNbIter = 1;
    for (int ii=0;ii<l_sNbIter;ii++)
    {
#endif
    // Compute the range [min,max] of the 2D image. This range will be divided into D_QUANTIZATION_RANGE_2D levels
    GetMinMaxVal_RGB(&l_tImage_min, &l_tImage_max, m_DataFusion.m_pt_I_R_DS, m_DataFusion.m_pt_I_G_DS, m_DataFusion.m_pt_I_B_DS, m_DataFusion.m_h*m_DataFusion.m_w);
    l_tImg_Range = l_tImage_max-l_tImage_min;
    s32 *l_pt_Jk_R[2];
    s32 *l_pl_Jk_R;
    s32 *l_pt_Jk_G[2];
    s32 *l_pl_Jk_G;
    s32 *l_pt_Jk_B[2];
    s32 *l_pl_Jk_B;
    s32 *l_pt_Jk_GRAY[2];
    s32 *l_pl_Jk_GRAY;
    s16 jk_0 = 0;
    s16 jk_1 = 1;
    s32 y,x;
    u08  l_ucShiftUnit = (1<<m_DataFusion.m_nr_shift);  // 1..16
#ifdef	MLF_BackGround_Counter
    s32 *l_pt_Jk_R_BG[2];
    s32 *l_pl_Jk_R_BG;
    l_pt_Jk_R_BG[0] = m_DataFusion.m_pt_Jk_level0_R_BG;
    l_pt_Jk_R_BG[1] = m_DataFusion.m_pt_Jk_level1_R_BG;
    s32 *l_pt_Jk_G_BG[2];
    s32 *l_pl_Jk_G_BG;
    l_pt_Jk_G_BG[0] = m_DataFusion.m_pt_Jk_level0_G_BG;
    l_pt_Jk_G_BG[1] = m_DataFusion.m_pt_Jk_level1_G_BG;
    s32 *l_pt_Jk_B_BG[2];
    s32 *l_pl_Jk_B_BG;
    l_pt_Jk_B_BG[0] = m_DataFusion.m_pt_Jk_level0_B_BG;
    l_pt_Jk_B_BG[1] = m_DataFusion.m_pt_Jk_level1_B_BG;
    s32 *l_pt_Jk_GRAY_BG[2];
    s32 *l_pl_Jk_GRAY_BG;
    l_pt_Jk_GRAY_BG[0] = m_DataFusion.m_pt_Jk_level0_GRAY_BG;
    l_pt_Jk_GRAY_BG[1] = m_DataFusion.m_pt_Jk_level1_GRAY_BG;
#endif
    l_pt_Jk_R[0] = m_DataFusion.m_pt_Jk_level0_R;
    l_pt_Jk_R[1] = m_DataFusion.m_pt_Jk_level1_R;
    l_pt_Jk_G[0] = m_DataFusion.m_pt_Jk_level0_G;
    l_pt_Jk_G[1] = m_DataFusion.m_pt_Jk_level1_G;
    l_pt_Jk_B[0] = m_DataFusion.m_pt_Jk_level0_B;
    l_pt_Jk_B[1] = m_DataFusion.m_pt_Jk_level1_B;
    l_pt_Jk_GRAY[0] = m_DataFusion.m_pt_Jk_level0_GRAY;
    l_pt_Jk_GRAY[1] = m_DataFusion.m_pt_Jk_level1_GRAY;

    // For all levels
    for (s16 i=0;i<D_QUANTIZATION_RANGE_2D;i++)
    {
        if (i==0)
        {
            l_pl_Jk_R = l_pt_Jk_R[jk_0];
            l_pl_Jk_G = l_pt_Jk_G[jk_0];
            l_pl_Jk_B = l_pt_Jk_B[jk_0];
            l_pl_Jk_GRAY = l_pt_Jk_GRAY[jk_0];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_R_BG = l_pt_Jk_R_BG[jk_0];
            l_pl_Jk_G_BG = l_pt_Jk_G_BG[jk_0];
            l_pl_Jk_B_BG = l_pt_Jk_B_BG[jk_0];
            l_pl_Jk_GRAY_BG = l_pt_Jk_GRAY_BG[jk_0];
#endif
        }
        else
        {
            l_pl_Jk_R = l_pt_Jk_R[jk_1];
            l_pl_Jk_G = l_pt_Jk_G[jk_1];
            l_pl_Jk_B = l_pt_Jk_B[jk_1];
            l_pl_Jk_GRAY = l_pt_Jk_GRAY[jk_1];
#ifdef	MLF_BackGround_Counter
            l_pl_Jk_R_BG = l_pt_Jk_R_BG[jk_1];
            l_pl_Jk_G_BG = l_pt_Jk_G_BG[jk_1];
            l_pl_Jk_B_BG = l_pt_Jk_B_BG[jk_1];
            l_pl_Jk_GRAY_BG = l_pt_Jk_GRAY_BG[jk_1];
#endif
        }
        {
            s32* l_pl_Wk_R = m_DataFusion.m_pt_Wk_R;
            s32* l_pl_Wk_G = m_DataFusion.m_pt_Wk_G;
            s32* l_pl_Wk_B = m_DataFusion.m_pt_Wk_B;
            s32* l_pl_Wk_GRAY = m_DataFusion.m_pt_Wk_GRAY;
            s32* l_pl_Jk_R_FG = l_pl_Jk_R;
            s32* l_pl_Jk_G_FG = l_pl_Jk_G;
            s32* l_pl_Jk_B_FG = l_pl_Jk_B;
            s32* l_pl_Jk_GRAY_FG = l_pl_Jk_GRAY;
#ifdef	MLF_BackGround_Counter
            s32* l_pl_Jk_R_BG2 = l_pl_Jk_R_BG;
            s32* l_pl_Jk_G_BG2 = l_pl_Jk_G_BG;
            s32* l_pl_Jk_B_BG2 = l_pl_Jk_B_BG;
            s32* l_pl_Jk_GRAY_BG2 = l_pl_Jk_GRAY_BG;
#endif
            // RGB to Grayscale
            s32 l_lIndexIni;
            if (i== 0)
                l_lIndexIni = l_tImage_min;
            else
            {
                if (i == l_ucQuantizationLevels)
                    l_lIndexIni = l_tImage_max;
                else
                    l_lIndexIni = ((s32)l_tImage_min +(i*(s32)(l_tImg_Range))/l_ucQuantizationLevels);
            }
            l_pt_D_DS = m_DataFusion.m_pt_D_DS;
            l_pt_I_R_DS = m_DataFusion.m_pt_I_R_DS;
            l_pt_I_G_DS = m_DataFusion.m_pt_I_G_DS;
            l_pt_I_B_DS = m_DataFusion.m_pt_I_B_DS;
            l_pt_I_GRAY_DS = m_DataFusion.m_pt_I_GRAY_DS;
            l_pt_Q_DS = m_DataFusion.m_pt_Q_DS;

            for (y=0;y<m_DataFusion.m_h;y++)
            {
                for (x=0;x<m_DataFusion.m_w;x++)
                {
                    s32 index_R=(l_lIndexIni-(s32)(*l_pt_I_R_DS++));
                    s32 index_G=(l_lIndexIni-(s32)(*l_pt_I_G_DS++));
                    s32 index_B=(l_lIndexIni-(s32)(*l_pt_I_B_DS++));
                    s32 index_GRAY=(l_lIndexIni-(s32)(*l_pt_I_GRAY_DS++));
                    if(index_R < 0)
                        index_R = -index_R;
                    if(index_G < 0)
                        index_G = -index_G;
                    if(index_B < 0)
                        index_B = -index_B;
                    if(index_GRAY < 0)
                        index_GRAY = -index_GRAY;
#ifdef	MLF_BackGround_Counter
                    if (*l_pt_D_DS >= MLF_APP_MAX_DISTANCE)
                    {
                        *l_pl_Wk_R++ = 0x0;
                        *l_pl_Jk_R_FG++ = 0x0;
                        *l_pl_Jk_R_BG2++ = (s32)(m_DataFusion.m_ucaTable[index_R]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Wk_G++ = 0x0;
                        *l_pl_Jk_G_FG++ = 0x0;
                        *l_pl_Jk_G_BG2++ = (s32)(m_DataFusion.m_ucaTable[index_G]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Wk_B++ = 0x0;
                        *l_pl_Jk_B_FG++ = 0x0;
                        *l_pl_Jk_B_BG2++ = (s32)(m_DataFusion.m_ucaTable[index_B]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Wk_GRAY++ = 0x0;
                        *l_pl_Jk_GRAY_FG++ = 0x0;
                        *l_pl_Jk_GRAY_BG2++ = (s32)(m_DataFusion.m_ucaTable[index_GRAY]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    }
                    else
                    {
                        *l_pl_Jk_R_BG2++ = 0x0;
                        *l_pl_Wk_R = (s32)(m_DataFusion.m_ucaTable[index_R]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_R_FG++ = (*l_pl_Wk_R)*(*l_pt_D_DS);
                        l_pl_Wk_R++;
                        *l_pl_Jk_G_BG2++ = 0x0;
                        *l_pl_Wk_G = (s32)(m_DataFusion.m_ucaTable[index_G]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_G_FG++ = (*l_pl_Wk_G)*(*l_pt_D_DS);
                        l_pl_Wk_G++;
                        *l_pl_Jk_B_BG2++ = 0x0;
                        *l_pl_Wk_B = (s32)(m_DataFusion.m_ucaTable[index_B]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_B_FG++ = (*l_pl_Wk_B)*(*l_pt_D_DS);
                        l_pl_Wk_B++;
                        *l_pl_Jk_GRAY_BG2++ = 0x0;
                        *l_pl_Wk_GRAY = (s32)(m_DataFusion.m_ucaTable[index_GRAY]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                        *l_pl_Jk_GRAY_FG++ = (*l_pl_Wk_GRAY)*(*l_pt_D_DS);
                        l_pl_Wk_GRAY++;
                    }
#else
                    *l_pl_Wk_R = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_R_FG++ = (*l_pl_Wk_R)*(*l_pt_D_DS);
                    l_pl_Wk_R++;
                    *l_pl_Wk_G = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_G_FG++ = (*l_pl_Wk_G)*(*l_pt_D_DS);
                    l_pl_Wk_G++;
                    *l_pl_Wk_B = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_B_FG++ = (*l_pl_Wk_B)*(*l_pt_D_DS);
                    l_pl_Wk_B++;
                    *l_pl_Wk_GRAY = (s32)(m_DataFusion.m_ucaTable[index]*(*l_pt_Q_DS)); // 255 Increase Accuracy
                    *l_pl_Jk_GRAY_FG++ = (*l_pl_Wk_GRAY)*(*l_pt_D_DS);
                    l_pl_Wk_GRAY++;
#endif
                    l_pt_D_DS++;
                    l_pt_Q_DS++;
                }
            }
        }
        {
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_R,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_G,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_G,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_B,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_B,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_GRAY,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,m_DataFusion.m_pt_Wk_GRAY,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#ifdef	MLF_BackGround_Counter
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_R_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_G_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_B_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
            Gaussian_Recursive_Order0(&m_DataFusion.m_StGaussian,l_pl_Jk_GRAY_BG,m_DataFusion.m_pdaImgBox,m_DataFusion.m_h,m_DataFusion.m_w);
#endif
        }
        {
            const s32* l_pl_Wk_R = m_DataFusion.m_pt_Wk_R;
            const s32* l_pl_Wk_G = m_DataFusion.m_pt_Wk_G;
            const s32* l_pl_Wk_B = m_DataFusion.m_pt_Wk_B;
            const s32* l_pl_Wk_GRAY = m_DataFusion.m_pt_Wk_GRAY;

            s32* l_pus_J_R = l_pl_Jk_R;
            s32* l_pus_J_G = l_pl_Jk_G;
            s32* l_pus_J_B = l_pl_Jk_B;
            s32* l_pus_J_GRAY = l_pl_Jk_GRAY;
#ifdef	MLF_BackGround_Counter
            s32* l_pus_J_R_BG = l_pl_Jk_R_BG;
            s32* l_pus_J_G_BG = l_pl_Jk_G_BG;
            s32* l_pus_J_B_BG = l_pl_Jk_B_BG;
            s32* l_pus_J_GRAY_BG = l_pl_Jk_GRAY_BG;
#endif
            y = (s32)m_DataFusion.m_h*(s32)m_DataFusion.m_w;
            while (y--)
            {
                {
#ifdef	MLF_BackGround_Counter
                    if (*l_pl_Wk_R <= *l_pus_J_R_BG) // FGa, in the case where the foreground weight is 0, we set the background mask to 1
                        *l_pus_J_R_BG = 255;
                    else
                        *l_pus_J_R_BG = (255*(*l_pus_J_R_BG))/(2*(*l_pl_Wk_R));
                    if (*l_pl_Wk_G <= *l_pus_J_G_BG)
                        *l_pus_J_G_BG = 255;
                    else
                        *l_pus_J_G_BG = (255*(*l_pus_J_G_BG))/(2*(*l_pl_Wk_G));
                    if (*l_pl_Wk_B <= *l_pus_J_B_BG)
                        *l_pus_J_B_BG = 255;
                    else
                        *l_pus_J_B_BG = (255*(*l_pus_J_B_BG))/(2*(*l_pl_Wk_B));
                    if (*l_pl_Wk_GRAY <= *l_pus_J_GRAY_BG)
                        *l_pus_J_GRAY_BG = 255;
                    else
                        *l_pus_J_GRAY_BG = (255*(*l_pus_J_GRAY_BG))/(2*(*l_pl_Wk_GRAY));
                    l_pus_J_R_BG++;
                    l_pus_J_G_BG++;
                    l_pus_J_B_BG++;
                    l_pus_J_GRAY_BG++;
#endif
                    if (*l_pl_Wk_R)
                        *l_pus_J_R = ((*l_pus_J_R))/(*l_pl_Wk_R);
                    if (*l_pl_Wk_G)
                        *l_pus_J_G = ((*l_pus_J_G))/(*l_pl_Wk_G);
                    if (*l_pl_Wk_B)
                        *l_pus_J_B = ((*l_pus_J_B))/(*l_pl_Wk_B);
                    if (*l_pl_Wk_GRAY)
                        *l_pus_J_GRAY = ((*l_pus_J_GRAY))/(*l_pl_Wk_GRAY);
                    l_pl_Wk_R++;
                    l_pl_Wk_G++;
                    l_pl_Wk_B++;
                    l_pl_Wk_GRAY++;
                    l_pus_J_R++;
                    l_pus_J_G++;
                    l_pus_J_B++;
                    l_pus_J_GRAY++;
                }
            }
        }
        // INTERPOLATION: Check which pixels in the high-resolution 2D image have an intensity value that belongs to the range (levels 0-1) that we've processed
        //                in order to update the final pixels
        if (i>0)
        {
            const s32 l_lXMax = (s32)((m_DataFusion.m_w-1)*l_ucShiftUnit);
            const s32 l_lYMax = (s32)((m_DataFusion.m_h-1)*l_ucShiftUnit)-1;
#ifdef	MLF_BackGround_Counter
            const s32* p_plImgJk0_BG;
            const s32* p_plImgJk0_R_BG = l_pt_Jk_R_BG[jk_0];
            const s32* p_plImgJk0_G_BG = l_pt_Jk_G_BG[jk_0];
            const s32* p_plImgJk0_B_BG = l_pt_Jk_B_BG[jk_0];
            const s32* p_plImgJk0_GRAY_BG = l_pt_Jk_GRAY_BG[jk_0];
            const s32* p_plImgJk0_d_BG;
            const s32* p_plImgJk0_d_R_BG = p_plImgJk0_R_BG+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_G_BG = p_plImgJk0_G_BG+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_B_BG = p_plImgJk0_B_BG+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_GRAY_BG = p_plImgJk0_GRAY_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_BG;
            const s32* p_plImgJk1_R_BG = l_pt_Jk_R_BG[jk_1];
            const s32* p_plImgJk1_G_BG = l_pt_Jk_G_BG[jk_1];
            const s32* p_plImgJk1_B_BG = l_pt_Jk_B_BG[jk_1];
            const s32* p_plImgJk1_GRAY_BG = l_pt_Jk_GRAY_BG[jk_1];
            const s32* p_plImgJk1_d_BG;
            const s32* p_plImgJk1_d_R_BG = p_plImgJk1_R_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_G_BG = p_plImgJk1_G_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_B_BG = p_plImgJk1_B_BG+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_GRAY_BG = p_plImgJk1_GRAY_BG+m_DataFusion.m_w;
#endif
            const s32* p_plImgJk0;
            const s32* p_plImgJk0_R = l_pt_Jk_R[jk_0];
            const s32* p_plImgJk0_G = l_pt_Jk_G[jk_0];
            const s32* p_plImgJk0_B = l_pt_Jk_B[jk_0];
            const s32* p_plImgJk0_GRAY = l_pt_Jk_GRAY[jk_0];
            const s32* p_plImgJk0_d;
            const s32* p_plImgJk0_d_R = p_plImgJk0_R+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_G = p_plImgJk0_G+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_B = p_plImgJk0_B+m_DataFusion.m_w;
            const s32* p_plImgJk0_d_GRAY = p_plImgJk0_GRAY+m_DataFusion.m_w;
            const s32* p_plImgJk1;
            const s32* p_plImgJk1_R = l_pt_Jk_R[jk_1];
            const s32* p_plImgJk1_G = l_pt_Jk_G[jk_1];
            const s32* p_plImgJk1_B = l_pt_Jk_B[jk_1];
            const s32* p_plImgJk1_GRAY = l_pt_Jk_GRAY[jk_1];
            const s32* p_plImgJk1_d;
            const s32* p_plImgJk1_d_R = p_plImgJk1_R+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_G = p_plImgJk1_G+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_B = p_plImgJk1_B+m_DataFusion.m_w;
            const s32* p_plImgJk1_d_GRAY = p_plImgJk1_GRAY+m_DataFusion.m_w;
            const s16* l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[0][0][0];
            const s32 l_lImgRangeXSihiftUnitSq = (s32)l_tImg_Range* (s32)l_ucShiftUnit*(s32)l_ucShiftUnit;
            const u08 *l_pt_I_R = m_DataFusion.m_pt_I_R;
            const u08 *l_pt_I_G = m_DataFusion.m_pt_I_G;
            const u08 *l_pt_I_B = m_DataFusion.m_pt_I_B;
            const u08 *l_pt_I_GRAY = m_DataFusion.m_pt_I_GRAY;
            const u16 *l_pt_D = m_DataFusion.m_pt_D;
            const u08 *l_pt_BetaChannel = m_DataFusion.m_pt_BetaChannel;
            const u08 *l_pt_BetaValue = m_DataFusion.m_pt_BetaValue;
            const u16 l_tKf1Level = (((i-1)*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            const u16 l_tKf2Level = (((i  )*(s16)l_tImg_Range)/(s16)l_ucQuantizationLevels)+l_tImage_min;
            s16 l_sY0 = 0;
            s16 l_sdY = 0;
            u16* l_pt_J = m_DataFusion.m_pt_J;

            for (y=0;y<m_DataFusion.m_h_original;y++)
            {
                {
                    s16 l_sdX = 0;
                    u08 l_ucVal;

                    for (x=0;x<l_lXMax;x++)
                    {
                        if (*l_pt_BetaChannel == D_CHANNEL_R) // R-Channel
                        {
                            l_ucVal = *l_pt_I_R;
                            p_plImgJk0_BG = p_plImgJk0_R_BG;
                            p_plImgJk0_d_BG = p_plImgJk0_d_R_BG;
                            p_plImgJk1_BG = p_plImgJk1_R_BG;
                            p_plImgJk1_d_BG = p_plImgJk1_d_R_BG;
                            p_plImgJk0 = p_plImgJk0_R;
                            p_plImgJk0_d = p_plImgJk0_d_R;
                            p_plImgJk1 = p_plImgJk1_R;
                            p_plImgJk1_d = p_plImgJk1_d_R;

                        }
                        else if (*l_pt_BetaChannel == D_CHANNEL_G) // G-Channel
                        {
                            l_ucVal = *l_pt_I_G;
                            p_plImgJk0_BG = p_plImgJk0_G_BG;
                            p_plImgJk0_d_BG = p_plImgJk0_d_G_BG;
                            p_plImgJk1_BG = p_plImgJk1_G_BG;
                            p_plImgJk1_d_BG = p_plImgJk1_d_G_BG;
                            p_plImgJk0 = p_plImgJk0_G;
                            p_plImgJk0_d = p_plImgJk0_d_G;
                            p_plImgJk1 = p_plImgJk1_G;
                            p_plImgJk1_d = p_plImgJk1_d_G;
                        }
                        else if (*l_pt_BetaChannel == D_CHANNEL_B) // B-Channel
                        {
                            l_ucVal = *l_pt_I_B;
                            p_plImgJk0_BG = p_plImgJk0_B_BG;
                            p_plImgJk0_d_BG = p_plImgJk0_d_B_BG;
                            p_plImgJk1_BG = p_plImgJk1_B_BG;
                            p_plImgJk1_d_BG = p_plImgJk1_d_B_BG;
                            p_plImgJk0 = p_plImgJk0_B;
                            p_plImgJk0_d = p_plImgJk0_d_B;
                            p_plImgJk1 = p_plImgJk1_B;
                            p_plImgJk1_d = p_plImgJk1_d_B;
                        }
                        else // GRAY-Channel
                        {
                            l_ucVal = *l_pt_I_GRAY;
                            p_plImgJk0_BG = p_plImgJk0_GRAY_BG;
                            p_plImgJk0_d_BG = p_plImgJk0_d_GRAY_BG;
                            p_plImgJk1_BG = p_plImgJk1_GRAY_BG;
                            p_plImgJk1_d_BG = p_plImgJk1_d_GRAY_BG;
                            p_plImgJk0 = p_plImgJk0_GRAY;
                            p_plImgJk0_d = p_plImgJk0_d_GRAY;
                            p_plImgJk1 = p_plImgJk1_GRAY;
                            p_plImgJk1_d = p_plImgJk1_d_GRAY;
                        }
                        if (l_ucVal<l_tKf2Level)
                        {
                           if (l_ucVal>=l_tKf1Level)
                           {
                                // qx_linear_interpolate_xy
                                s32 l_lRetInterpolate2;
                                s32 l_lRetInterpolate;
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d_BG+1));
                                l_lRetInterpolate *= alpha;
                                l_psaDxx -= 3; // reset the position
                                l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate2 *= alpha_Comp;
                                l_psaDxx -= 3; // reset the position
                                l_BakGroundWeight =(u16)((l_lRetInterpolate+l_lRetInterpolate2) /(l_lImgRangeXSihiftUnitSq));
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk0));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk0_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk0_d+1));
                                    l_lRetInterpolate *= alpha;

                                    l_psaDxx -= 3;// reset the position
                                    l_lRetInterpolate2  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate2 += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate2 += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate2 *= alpha_Comp;

                                    if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                        *l_pt_J =((u16)(255-(*l_pt_BetaValue))*((l_lRetInterpolate+l_lRetInterpolate2)/(l_lImgRangeXSihiftUnitSq)) + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                    else
                                        *l_pt_J =(u16)((l_lRetInterpolate+l_lRetInterpolate2)/(l_lImgRangeXSihiftUnitSq));
                                }
                            }
                        }
                        else
                        {
                            if (i==l_ucQuantizationLevels)
                            {
                                const s16* l_psaDxx = l_spaDxx_Ylevel+l_sdX*4;
                                s32 l_lRetInterpolate ;
#ifdef	MLF_BackGround_Counter
                                s32 l_BakGroundWeight;
                                l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1_BG));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_BG+1));
                                l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d_BG));
                                l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d_BG+1));
                                l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                l_psaDxx -= 3;// reset the position
                                l_BakGroundWeight = l_lRetInterpolate;
                                if (l_BakGroundWeight > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    l_lRetInterpolate  = (*l_psaDxx++)    *(*(p_plImgJk1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1+1));
                                    l_lRetInterpolate += (*l_psaDxx++)    *(*(p_plImgJk1_d));
                                    l_lRetInterpolate += (*l_psaDxx)      *(*(p_plImgJk1_d+1));
                                    l_lRetInterpolate += (l_ucShiftUnit*l_ucShiftUnit)/2/*0.5f*/;
                                    l_lRetInterpolate /= (l_ucShiftUnit*l_ucShiftUnit);
                                    if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                        *l_pt_J =((u16)(255-(*l_pt_BetaValue))*l_lRetInterpolate + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                    else
                                        *l_pt_J =(u16)(l_lRetInterpolate);
                                }
                            }
                        }
                        l_pt_J++;
                        l_pt_I_R++;
                        l_pt_I_B++;
                        l_pt_I_G++;
                        l_pt_I_GRAY++;
                        l_pt_D++;
                        l_pt_BetaChannel++;
                        l_pt_BetaValue++;
                        {
                            l_sdX++;
                            if (l_sdX == l_ucShiftUnit)
                            {
                                l_sdX= 0;
#ifdef	MLF_BackGround_Counter
                                p_plImgJk0_R_BG++;
                                p_plImgJk0_G_BG++;
                                p_plImgJk0_B_BG++;
                                p_plImgJk0_GRAY_BG++;
                                p_plImgJk0_d_R_BG++;
                                p_plImgJk0_d_G_BG++;
                                p_plImgJk0_d_B_BG++;
                                p_plImgJk0_d_GRAY_BG++;
                                p_plImgJk1_R_BG++;
                                p_plImgJk1_G_BG++;
                                p_plImgJk1_B_BG++;
                                p_plImgJk1_GRAY_BG++;
                                p_plImgJk1_d_R_BG++;
                                p_plImgJk1_d_G_BG++;
                                p_plImgJk1_d_B_BG++;
                                p_plImgJk1_d_GRAY_BG++;
#endif

                                p_plImgJk0_R++;
                                p_plImgJk0_G++;
                                p_plImgJk0_B++;
                                p_plImgJk0_GRAY++;
                                p_plImgJk0_d_R++;
                                p_plImgJk0_d_G++;
                                p_plImgJk0_d_B++;
                                p_plImgJk0_d_GRAY++;
                                p_plImgJk1_R++;
                                p_plImgJk1_G++;
                                p_plImgJk1_B++;
                                p_plImgJk1_GRAY++;
                                p_plImgJk1_d_R++;
                                p_plImgJk1_d_G++;
                                p_plImgJk1_d_B++;
                                p_plImgJk1_d_GRAY++;
                            }
                        }
                    }
                    // X Max to end
                    {
                        s32 l_dLIc;
                        s32 l_dLIa;
                        s32 l_dLIb;
                        const s16* l_psaDxx = l_spaDxx_Ylevel;
                        /*qx_linear_interpolate_xy */
#ifdef	MLF_BackGround_Counter
                        s32 l_dLIc_BG;
                        s32 l_dLIa_BG;
                        s32 l_dLIb_BG;
                        l_dLIa_BG = ((*l_psaDxx)  *(*(p_plImgJk0_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d_BG)));
                        l_dLIb_BG = (  (*l_psaDxx)*(*(p_plImgJk1_BG))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d_BG)));
                        l_dLIc_BG = l_dLIb_BG+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
#endif
                        l_dLIa = ((*l_psaDxx)  *(*(p_plImgJk0))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk0_d)));
                        l_dLIb = (  (*l_psaDxx)*(*(p_plImgJk1))+
                                 (*(l_psaDxx+2))*(*(p_plImgJk1_d)));
                        l_dLIc = l_dLIb+(l_ucShiftUnit*l_ucShiftUnit)/2/*+*0.5f*/;
                        for (x = l_lXMax;x<m_DataFusion.m_w_original;x++)
                        {
                            if (l_ucVal<l_tKf2Level)
                            {
                              if (l_ucVal>=l_tKf1Level)
                              {
                                const s16 alpha= (l_tKf2Level - l_ucVal)* l_ucQuantizationLevels; // 0 <=alpha<=l_tImg_Range
                                const s16 alpha_Comp = (l_tImg_Range-alpha);

#ifdef	MLF_BackGround_Counter
                                *l_pt_J =(u16)(((alpha*l_dLIa_BG+(alpha_Comp)* l_dLIb_BG))/(l_lImgRangeXSihiftUnitSq));
                                if (*l_pt_J > D_BG_THRESHOLD)
                                    *l_pt_J = MLF_APP_MAX_DISTANCE;
                                else
#endif
                                {
                                    if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                        *l_pt_J =((u16)(255-(*l_pt_BetaValue))* (((alpha*l_dLIa+(alpha_Comp)*l_dLIb))/(l_lImgRangeXSihiftUnitSq)) + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                    else
                                        *l_pt_J =(u16)(((alpha*l_dLIa+(alpha_Comp)* l_dLIb))/(l_lImgRangeXSihiftUnitSq));
                                }
                              }
                            }
                            else
                            {
                                if (i==l_ucQuantizationLevels)
                                {
#ifdef	MLF_BackGround_Counter
                                    *l_pt_J =(u16)(l_dLIc_BG/(l_ucShiftUnit*l_ucShiftUnit));
                                    if (*l_pt_J > D_BG_THRESHOLD)
                                        *l_pt_J = MLF_APP_MAX_DISTANCE;
                                    else
#endif
                                        if (*l_pt_D != MLF_APP_MAX_DISTANCE)
                                            *l_pt_J =((u16)(255-(*l_pt_BetaValue))* (l_dLIc/(l_ucShiftUnit*l_ucShiftUnit)) + ((u16)*l_pt_BetaValue)*(*l_pt_D))/255;
                                        else
                                            *l_pt_J =(u16)(l_dLIc/(l_ucShiftUnit*l_ucShiftUnit));
                                }
                            }
                            l_pt_J++;
                            l_pt_I_R++;
                            l_pt_I_G++;
                            l_pt_I_B++;
                            l_pt_I_GRAY++;
                            l_pt_D++;
                            l_pt_BetaChannel++;
                            l_pt_BetaValue++;
                        }
                    }
                }
                if (y < l_lYMax)
                {
                    l_sdY++;
                    if (l_sdY == l_ucShiftUnit)
                    {
                        l_sdY= 0;
                        l_sY0++;
                    }
                    l_spaDxx_Ylevel = &m_DataFusion.m_saDxx[l_sdY][0][0];
                }
#ifdef	MLF_BackGround_Counter
                p_plImgJk0_R_BG   = l_pt_Jk_R_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_G_BG   = l_pt_Jk_G_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_B_BG   = l_pt_Jk_B_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_GRAY_BG   = l_pt_Jk_GRAY_BG[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_R_BG = p_plImgJk0_R_BG + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_G_BG = p_plImgJk0_G_BG + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_B_BG = p_plImgJk0_B_BG + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_GRAY_BG = p_plImgJk0_GRAY_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_R_BG   = l_pt_Jk_R_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_G_BG   = l_pt_Jk_G_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_B_BG   = l_pt_Jk_B_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_GRAY_BG   = l_pt_Jk_GRAY_BG[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_R_BG = p_plImgJk1_R_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_G_BG = p_plImgJk1_G_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_B_BG = p_plImgJk1_B_BG + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_GRAY_BG = p_plImgJk1_GRAY_BG + (s32)m_DataFusion.m_w;
#endif

                p_plImgJk0_R   = l_pt_Jk_R[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_G   = l_pt_Jk_G[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_B   = l_pt_Jk_B[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_GRAY   = l_pt_Jk_GRAY[jk_0] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk0_d_R = p_plImgJk0_R + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_G = p_plImgJk0_G + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_B = p_plImgJk0_B + (s32)m_DataFusion.m_w;
                p_plImgJk0_d_GRAY = p_plImgJk0_GRAY + (s32)m_DataFusion.m_w;
                p_plImgJk1_R   = l_pt_Jk_R[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_G   = l_pt_Jk_G[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_B   = l_pt_Jk_B[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_GRAY   = l_pt_Jk_GRAY[jk_1] + (s32)m_DataFusion.m_w * (s32)l_sY0;
                p_plImgJk1_d_R = p_plImgJk1_R + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_G = p_plImgJk1_G + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_B = p_plImgJk1_B + (s32)m_DataFusion.m_w;
                p_plImgJk1_d_GRAY = p_plImgJk1_GRAY + (s32)m_DataFusion.m_w;
            }
            jk_1=jk_0;
            jk_0=(jk_0+1)%2;
        }
    }
#ifdef MLF_ANALYSE_TIME_CONSUMPTION
    }
    std::clock_t c_end = std::clock();

    std::cout << "CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC
              << " ms" << std::endl;
    std::cout << "Done..." << std::endl;
#endif

    return 0;
}

/* Apply Joint Bilateral Upsampling Filter (JBU Filter)
 * J. Kopf, M. Cohen, D. Lischinski, M. Uttyttendaele. Joint bilateral upsampling. ACM Transactions on Graphics (TOG), 26(3), 2007
 */
short c_DataFusion::JBU_Filter_Kopf()
{
    u16* l_pt_J; // Enhanced depth map, J
    u16* l_pt_D_q; // q index within the downsampled depth map D_DS
    u08* l_pt_I_p; // p index in I
    u08* l_pt_I_q; // q index in I    
    f32 l_fRangeTerm;
    f32 l_fSigSpatial = m_DataFusion.m_Set.m_fSigma_spatial; // Sigma spatial
    f32 l_fSigRange = m_DataFusion.m_Set.m_fSigma_range; // Sigma range
    s16 l_sFilterSize = 12; // Filter support, neighborhood is 31x31
    f32 l_faSpatialTerm[25][25]; // Precomputed array for the spatial term of 31x31 values
    f32 l_fSum, l_fSumNorm, l_fAux;
    f32 l_fExpSpatial, l_fExpRange;
    f32 l_fDistSpatial, l_fDistRange;
    f32 l_fSigSpatialSquare, l_fSigRangeSquare;
    f32 l_fNormFactorSpatial, l_fNormFactorRange;

    l_fSigSpatialSquare = 2*l_fSigSpatial*l_fSigSpatial;
    l_fNormFactorSpatial = 1/(M_PI*l_fSigSpatialSquare);
    l_fSigRangeSquare = 2*l_fSigRange*l_fSigRange;
    l_fNormFactorRange = 1/(M_PI*l_fSigRangeSquare);

    std::cout << "Starting JBU Filter (Kopf et al.)" << std::endl;
    std::cout << "Precomputing values..." << std::endl;
    // Precompute spatial terms
    for (s16 jj=-l_sFilterSize;jj<=l_sFilterSize;jj++)
    {
        for (s16 ii=-l_sFilterSize;ii<=l_sFilterSize;ii++)
        {
            l_fDistSpatial = -(ii*ii+jj*jj);
            l_fExpSpatial = exp(l_fDistSpatial/l_fSigSpatialSquare);
            l_faSpatialTerm[jj+l_sFilterSize][ii+l_sFilterSize] = l_fNormFactorSpatial*l_fExpSpatial;
        }
    }
    std::cout << "Start filtering..." << std::endl;
    // The border is not filled
    for (s16 j=l_sFilterSize;j<m_DataFusion.m_h_original-l_sFilterSize;j++)
    {
        for (s16 i=l_sFilterSize;i<m_DataFusion.m_w_original-l_sFilterSize;i++)
        {
            l_pt_I_p = m_DataFusion.m_pt_I_GRAY + j*m_DataFusion.m_w_original + i; // Pointer to the guidance image, I_p
            l_pt_J = m_DataFusion.m_pt_J + j*m_DataFusion.m_w_original + i;        // Pointer to the enhanced depth map, J
            l_fSum = 0;     // Reset sum
            l_fSumNorm = 0; // Reset norm factor

            // Tackle the neighborhood, q
            for (s16 jj=-l_sFilterSize;jj<=l_sFilterSize;jj++)
            {
                for (s16 ii=-l_sFilterSize;ii<=l_sFilterSize;ii++)
                {
                    l_pt_I_q = m_DataFusion.m_pt_I_GRAY + (j+jj)*m_DataFusion.m_w_original + (i+ii); // Pointer to the guidance image, I_q
                    l_pt_D_q = m_DataFusion.m_pt_D + (j+jj)*m_DataFusion.m_w_original + (i+ii);   // Pointer to the downsampled depth map, D_DS_q

                    if (*l_pt_D_q == 0) // Bg. handling
                        continue;

                    l_fDistRange = -(std::pow((f32)((*l_pt_I_p)-(*l_pt_I_q)),2));
                    l_fExpRange = exp(l_fDistRange/l_fSigRangeSquare);
                    l_fRangeTerm = l_fNormFactorRange*l_fExpRange;
                    l_fAux = l_faSpatialTerm[jj+l_sFilterSize][ii+l_sFilterSize]*l_fRangeTerm;
                    l_fSumNorm += l_fAux;
                    l_fSum += l_fAux*(*l_pt_D_q);
                }
            }
            *l_pt_J = (u16)(l_fSum/l_fSumNorm); // Update value to J
        }
    }
    std::cout << "JBU Filter (Kopf et al.) done..." << std::endl;

    return 0;
}

/* Apply New Joint Bilateral Upsampling Filter
 * S. Y. Kim, J. H. Cho, A. Koschan, M. A. Abidi. Spatial and Temporal Enhancement of Depth Images Captured by a Time-of-Flight Depth Sensor
 */
short c_DataFusion::NJBU_Filter_Kim()
{
    u16* l_pt_J; // Enhanced depth map, J
    u16* l_pt_D_p; // p index within the downsampled depth map D
    u16* l_pt_D_q; // q index within the downsampled depth map D
    u08* l_pt_I_p; // p index in I
    u08* l_pt_I_q; // q index in I
    f32 l_fRangeTerm_I, l_fRangeTerm_D;
    f32 l_fSigSpatial = m_DataFusion.m_Set.m_fSigma_spatial; // Sigma spatial
    f32 l_fSigRange_I = m_DataFusion.m_Set.m_fSigma_range; // Sigma range
    f32 l_fSigRange_D = m_DataFusion.m_Set.m_fSigma_range; // Sigma range
    s16 l_sFilterSize = 12; // Filter support, neighborhood is 9x9 which is 3x3 in the low res image
    f32 l_faSpatialTerm[D_KERNEL_SIZE][D_KERNEL_SIZE]; // Precomputed array for the spatial term of 31x31 values
    f32 l_fSum, l_fSumNorm, l_fAlpha, l_fAux;
    f32 l_fDistSpatial, l_fDistRange_I, l_fDistRange_D;
    f32 l_fExpSpatial, l_fExpRange_I, l_fExpRange_D;
    f32 l_fSigSpatialSquare, l_fSigRangeSquare_I, l_fSigRangeSquare_D;
    f32 l_fNormFactorSpatial, l_fNormFactorRange_I, l_fNormFactorRange_D;

    l_fSigSpatialSquare = 2*l_fSigSpatial*l_fSigSpatial;
    l_fNormFactorSpatial = 1/(M_PI*l_fSigSpatialSquare);
    l_fSigRangeSquare_I = 2*l_fSigRange_I*l_fSigRange_I;
    l_fNormFactorRange_I = 1/(M_PI*l_fSigRangeSquare_I);
    l_fSigRangeSquare_D = 2*l_fSigRange_D*l_fSigRange_D;
    l_fNormFactorRange_D = 1/(M_PI*l_fSigRangeSquare_D);

    std::cout << "Starting NJBU Filter (Kim et al.)" << std::endl;
    std::cout << "Precomputing values..." << std::endl;
    // Precompute spatial terms
    for (s16 jj=-l_sFilterSize;jj<=l_sFilterSize;jj++)
    {
        for (s16 ii=-l_sFilterSize;ii<=l_sFilterSize;ii++)
        {
            l_fDistSpatial = -(ii*ii+jj*jj);
            l_fExpSpatial = exp(l_fDistSpatial/l_fSigSpatialSquare);
            l_faSpatialTerm[jj+l_sFilterSize][ii+l_sFilterSize] = l_fNormFactorSpatial*l_fExpSpatial;
        }
    }

    std::cout << "Start filtering..." << std::endl;
    // The border is not filled
    for (s16 j=l_sFilterSize;j<m_DataFusion.m_h_original-l_sFilterSize;j++)
    {
        for (s16 i=l_sFilterSize;i<m_DataFusion.m_w_original-l_sFilterSize;i++)
        {
            l_pt_I_p = m_DataFusion.m_pt_I_GRAY + j*m_DataFusion.m_w_original + i; // Pointer to the guidance image, I_p
            l_pt_D_p = m_DataFusion.m_pt_D + j*m_DataFusion.m_w_original + i;   // Pointer to the downsampled depth map, D_DS_p
            l_pt_J = m_DataFusion.m_pt_J + j*m_DataFusion.m_w_original + i;        // Pointer to the enhanced depth map, J
            l_fSum = 0;     // Reset sum
            l_fSumNorm = 0; // Reset norm factor

            // Tackle the neighborhood, q
            for (s16 jj=-l_sFilterSize;jj<=l_sFilterSize;jj++)
            {
                for (s16 ii=-l_sFilterSize;ii<=l_sFilterSize;ii++)
                {
                    l_pt_I_q = m_DataFusion.m_pt_I_GRAY + (j+jj)*m_DataFusion.m_w_original + (i+ii); // Pointer to the guidance image, I_q
                    l_pt_D_q = m_DataFusion.m_pt_D + (j+jj)*m_DataFusion.m_w_original + (i+ii);   // Pointer to the downsampled depth map, D_DS_q

                    if (*l_pt_D_q == 0) // Bg. handling
                        continue;

                    l_fDistRange_I = -(std::pow((f32)((*l_pt_I_p)-(*l_pt_I_q)),2));
                    l_fExpRange_I = exp(l_fDistRange_I/l_fSigRangeSquare_I);
                    l_fDistRange_D = -(std::pow((f32)((*l_pt_D_p)-(*l_pt_D_q)),2));
                    l_fExpRange_D = exp(l_fDistRange_D/l_fSigRangeSquare_D);
                    l_fRangeTerm_I = l_fNormFactorRange_I*l_fExpRange_I;
                    l_fRangeTerm_D = l_fNormFactorRange_D*l_fExpRange_D;
                    l_fAux = l_faSpatialTerm[jj+l_sFilterSize][ii+l_sFilterSize]*l_fRangeTerm_I*l_fRangeTerm_D;
                    l_fSumNorm += l_fAux;
                    l_fSum += l_fAux*(*l_pt_D_q);
                }
            }
            *l_pt_J = (u16)(l_fSum/l_fSumNorm); // Update value to J
        }
    }
    std::cout << "NJBU Filter (Kim et al.) done..." << std::endl;

    return 0;
}

/* Apply Noise-Aware Filter (NAFDU Filter)
 * D. Chan, H. Buisman, C. Theobalt, S. Thrun. A Noise-Aware Filter for Real-Time Depth Upsampling. M2SFA2 2008.
 */
short c_DataFusion::NAFDU_Filter_Chan()
{
    u16* l_pt_J; // Enhanced depth map, J
    u16* l_pt_D_p; // p index within the downsampled depth map D
    u16* l_pt_D_q; // q index within the downsampled depth map D
    u08* l_pt_I_p; // p index in I
    u08* l_pt_I_q; // q index in I
    f32 l_fRangeTerm_I, l_fRangeTerm_D;
    f32 l_fSigSpatial = m_DataFusion.m_Set.m_fSigma_spatial; // Sigma spatial
    f32 l_fSigRange_I = m_DataFusion.m_Set.m_fSigma_range; // Sigma range
    f32 l_fSigRange_D = m_DataFusion.m_Set.m_fSigma_range; // Sigma range
    s16 l_sFilterSize = 12; // Filter support, neighborhood is 9x9 which is 3x3 in the low res image
    f32 l_faSpatialTerm[D_KERNEL_SIZE][D_KERNEL_SIZE]; // Precomputed array for the spatial term of 31x31 values
    f32 l_fSum, l_fSumNorm, l_fAlpha, l_fAux;
    f32 l_fDistSpatial, l_fDistRange_I, l_fDistRange_D;
    f32 l_fExpSpatial, l_fExpRange_I, l_fExpRange_D;
    f32 l_fEpsilon; // Controls how wide the transition are (in terms of min-max difference) between case I and case II is
    f32 l_fTau; // Controls at what min-max difference the blending interval shall be centered
    f32 l_fSigSpatialSquare, l_fSigRangeSquare_I, l_fSigRangeSquare_D;
    f32 l_fNormFactorSpatial, l_fNormFactorRange_I, l_fNormFactorRange_D;

    l_fSigSpatialSquare = 2*l_fSigSpatial*l_fSigSpatial;
    l_fNormFactorSpatial = 1/(M_PI*l_fSigSpatialSquare);
    l_fSigRangeSquare_I = 2*l_fSigRange_I*l_fSigRange_I;
    l_fNormFactorRange_I = 1/(M_PI*l_fSigRangeSquare_I);
    l_fSigRangeSquare_D = 2*l_fSigRange_D*l_fSigRange_D;
    l_fNormFactorRange_D = 1/(M_PI*l_fSigRangeSquare_D);

    std::cout << "Starting NAFDU Filter (Chan et al.)" << std::endl;
    std::cout << "Precomputing values..." << std::endl;
    // Precompute spatial terms
    for (s16 jj=-l_sFilterSize;jj<=l_sFilterSize;jj++)
    {
        for (s16 ii=-l_sFilterSize;ii<=l_sFilterSize;ii++)
        {
            l_fDistSpatial = -(ii*ii+jj*jj);
            l_fExpSpatial = exp(l_fDistSpatial/l_fSigSpatialSquare);
            l_faSpatialTerm[jj+l_sFilterSize][ii+l_sFilterSize] = l_fNormFactorSpatial*l_fExpSpatial;
        }
    }

    // Dilate, erode and substract the input image to have the max-min (inc of omega)
    cv::Mat l_cvElement = cv::getStructuringElement(0, cv::Size(D_KERNEL_SIZE, D_KERNEL_SIZE));
    cv::Mat l_cvD(cv::Size(m_DataFusion.m_w_original, m_DataFusion.m_h_original), CV_16UC1, (unsigned char*)m_DataFusion.m_pt_D, cv::Mat::AUTO_STEP);
    cv::Mat l_cvD_Blurred;
    // Low pas filter (Gaussian blur) to detect only noise affected regions
    cv::blur(l_cvD, l_cvD_Blurred, cv::Size(3, 3));
    cv::Mat l_cvEroded; // Min value
    cv::erode(l_cvD_Blurred, l_cvEroded, l_cvElement);
    cv::Mat l_cvDilated; // Max value
    cv::dilate(l_cvD_Blurred, l_cvDilated, l_cvElement);
    cv::Mat l_cvDiffDilateErode;
    l_cvDiffDilateErode = l_cvDilated-l_cvEroded;
    u16* l_pt_DiffOmega = (u16*)l_cvDiffDilateErode.data; // q index within the downsampled depth map D
    double l_dMinDiffOmega, l_dMaxDiffOmega;
    cv::minMaxLoc(l_cvDiffDilateErode, &l_dMinDiffOmega, &l_dMaxDiffOmega); // Find minimum and maximum inc omega
    l_fTau = ((f32)l_dMaxDiffOmega-(f32)l_dMinDiffOmega)/2.0;
    l_fEpsilon = 1/(m_DataFusion.m_Set.m_fSigma_spatial*m_DataFusion.m_Set.m_fSigma_spatial);
/*
    //Start debug
    // i=285, j=83
    s16 l_s_i=382, l_s_j = 7;
    s32  l_lIndex = l_s_j*m_DataFusion.m_w_original + l_s_i;
    u16* l_ptD_data = (u16*)l_cvD.data + l_lIndex;
    u16* l_ptD_blurred_data = (u16*)l_cvD_Blurred.data + l_lIndex;
    u16* l_ptMax_D_data = (u16*)l_cvDilated.data + l_lIndex;
    u16* l_ptMin_D_data = (u16*)l_cvEroded.data + l_lIndex;
    l_pt_DiffOmega += l_lIndex;
    s32  l_Max=0, l_Min=9999;
    u16* l_ptD_data_window;

    for (s16 jj=-l_sFilterSize;jj<=l_sFilterSize;jj++)
    {
        for (s16 ii=-l_sFilterSize;ii<=l_sFilterSize;ii++)
        {
            l_ptD_data_window = (u16*)l_cvD.data + (l_s_j+jj)*m_DataFusion.m_w_original + (l_s_i+ii);
            if (l_Max < *l_ptD_data_window)
                l_Max = *l_ptD_data_window;
            if (l_Min > *l_ptD_data_window)
                l_Min = *l_ptD_data_window;
        }
    }
    // Cout
    std::cout << "D(250,250): " << *l_ptD_data << ", denoised: " << *l_ptD_blurred_data << std::endl;
    std::cout << "Max value: " << l_Max << ", Dilated: " << *l_ptMax_D_data << std::endl;
    std::cout << "Min value: " << l_Min << ", Eroded: " << *l_ptMin_D_data << std::endl;
    std::cout << "Diff value: " << (l_Max-l_Min) << ", l_pt_DiffOmega: " << *l_pt_DiffOmega << std::endl;

    return 0;
    //End
*/
    std::cout << "Start filtering..." << std::endl;
    // The border is not filled
    for (s16 j=l_sFilterSize;j<m_DataFusion.m_h_original-l_sFilterSize;j++)
    {
        for (s16 i=l_sFilterSize;i<m_DataFusion.m_w_original-l_sFilterSize;i++)
        {
            l_pt_I_p = m_DataFusion.m_pt_I_GRAY + j*m_DataFusion.m_w_original + i; // Pointer to the guidance image, I_p
            l_pt_D_p = m_DataFusion.m_pt_D + j*m_DataFusion.m_w_original + i;   // Pointer to the downsampled depth map, D_DS_p
            l_pt_J = m_DataFusion.m_pt_J + j*m_DataFusion.m_w_original + i;        // Pointer to the enhanced depth map, J
            l_pt_DiffOmega = (u16*)l_cvDiffDilateErode.data + j*l_cvDiffDilateErode.cols + i;
            l_fSum = 0;     // Reset sum
            l_fSumNorm = 0; // Reset norm factor                    

            // Tackle the neighborhood, q
            for (s16 jj=-l_sFilterSize;jj<=l_sFilterSize;jj++)
            {
                for (s16 ii=-l_sFilterSize;ii<=l_sFilterSize;ii++)
                {
                    l_pt_I_q = m_DataFusion.m_pt_I_GRAY + (j+jj)*m_DataFusion.m_w_original + (i+ii); // Pointer to the guidance image, I_q
                    l_pt_D_q = m_DataFusion.m_pt_D + (j+jj)*m_DataFusion.m_w_original + (i+ii);   // Pointer to the downsampled depth map, D_DS_q

                    if (*l_pt_D_q == 0) // Bg. handling
                        continue;

                    l_fDistRange_I = -(std::pow((f32)((*l_pt_I_p)-(*l_pt_I_q)),2));
                    l_fExpRange_I = exp(l_fDistRange_I/l_fSigRangeSquare_I);
                    l_fDistRange_D = -(std::pow((f32)((*l_pt_D_p)-(*l_pt_D_q)),2));
                    l_fExpRange_D = exp(l_fDistRange_D/l_fSigRangeSquare_D);
                    l_fRangeTerm_I = l_fNormFactorRange_I*l_fExpRange_I;
                    l_fRangeTerm_D = l_fNormFactorRange_D*l_fExpRange_D;
                    l_fTau = (*l_pt_DiffOmega)/2;
                    l_fAlpha = 1/(1+exp(-l_fEpsilon*((*l_pt_DiffOmega)-l_fTau))); // Blending function, alpha
                    l_fAux = l_faSpatialTerm[jj+l_sFilterSize][ii+l_sFilterSize]*(l_fAlpha*l_fRangeTerm_I + (1-l_fAlpha)*l_fRangeTerm_D);
                    l_fSumNorm += l_fAux;
                    l_fSum += l_fAux*(*l_pt_D_q);
                }
            }
            *l_pt_J = (u16)(l_fSum/l_fSumNorm); // Update value to J
        }
    }
    std::cout << "NAFDU Filter (Chan et al.) done..." << std::endl;

    return 0;
}
