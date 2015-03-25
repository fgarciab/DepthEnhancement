/*
* Author: Frederic GARCIA BECERRO
* Email: frederic.garcia.becerro@gmail.com
* Website: http://www.frederic-garcia-becerro.com
*/

#ifndef C_DATAFUSION_H
#define C_DATAFUSION_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include "app_types.h"
#include <ctime>
// OpenCV
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"


#define D_QUANTIZATION_RANGE_3D		300 //160 /* It's related to the range data, sigmas, etc */
#define D_QUANTIZATION_RANGE_2D		50 //16
#define D_BG_THRESHOLD				128 // Background threshold // Originally 128, then shifted to 15, then back to 128
#define D_CHANNEL_GRAY				150 // Just used for display colours
#define D_CHANNEL_R					255
#define D_CHANNEL_G					128
#define D_CHANNEL_B					0
#define D_KERNEL_SIZE               25
#define MLF_APP_MAX_DISTANCE		20000   //20000	 // IEE ToF Max Dist: 20000
#define MLF_APP_2D_RANGE			256
#define MLF_APP_INTENSITY_RANGE		4096 /*define type size of image*/
#define MLF_BackGround_Counter      1

// FGa, from previous tool
#define D_Fusion_ScaleFactor        3   // Scale factor: 3x
#define D_Fusion_2D_3D_2DW			640 // By default we booked the maximum memory, considering VGA Res
#define D_Fusion_2D_3D_2DH			480
#define D_Fusion_2D_3D_2DSIZE		(D_Fusion_2D_3D_2DW*D_Fusion_2D_3D_2DH)
#define D_Fusion_2D_3D_2DSIZE_RGB	(D_Fusion_2D_3D_2DW*D_Fusion_2D_3D_2DH*3)

#define D_Fusion_2D_3D_Spatial_downsample_rate	0 // Default downsampling rate, can be [0,4] but I set it to 0 in order to book enough memory
#define D_Fusion_2D_3D_2DW_DS		(D_Fusion_2D_3D_2DW >> D_Fusion_2D_3D_Spatial_downsample_rate)
#define D_Fusion_2D_3D_2DH_DS		(D_Fusion_2D_3D_2DH >> D_Fusion_2D_3D_Spatial_downsample_rate)
#define D_Fusion_2D_3D_2DSIZE_DS	(D_Fusion_2D_3D_2DW_DS*D_Fusion_2D_3D_2DH_DS)

#define D_MEM_s32_DS		D_Fusion_2D_3D_2DSIZE_DS*sizeof(s32)
#define D_MEM_u16_DS		D_Fusion_2D_3D_2DSIZE_DS*sizeof(u16)
#define D_MEM_u08_DS		D_Fusion_2D_3D_2DSIZE_DS*sizeof(u08)
#define D_MEM_s32			D_Fusion_2D_3D_2DSIZE*sizeof(s32)
#define D_MEM_u16			D_Fusion_2D_3D_2DSIZE*sizeof(u16)
#define D_MEM_u08			D_Fusion_2D_3D_2DSIZE*sizeof(u08)
#define D_MEM_RGB_u08		D_Fusion_2D_3D_2DSIZE*3

#define D_Fusion_2D_3D_MEM_SIZE		(21*D_MEM_s32_DS + 2*D_MEM_u16_DS + 5*D_MEM_u08_DS + 2*D_MEM_u16 + 7*D_MEM_u08 + D_MEM_RGB_u08) // Totally required memory

// Filter type
typedef enum
{
    BF_FILTER       = 0, // Bilateral Filter
    JBU_FILTER      = 1, // JBU Filter
    PWAS_FILTER     = 2, // PWAS Filter
    PWAS_RGB_FILTER = 3, // PWAS Filter with I RGB
    UML_FILTER      = 4, // UML Filter, second PWAS neglected, values directly from high-res D   
    UML_RGB_FILTER  = 5, // UML Filter, second PWAS neglected, values directly from high-res D and I RGB
    RGBD_FILTER     = 6, // RGB-D Filter
    LITERATURE_FILTERS    = 7  // Literature Filters
} DataFusion_FilterType;

typedef struct
{
    s32 m_lMemSize;
    s32 m_lMemUsed;
    s08* m_pMemPointer;
}T_Memory_St;

typedef struct
{
    f32 m_fSigma_range;	// 0 ..1
    f32 m_fSigma_spatial;
    f32 m_fSigma_Q;
    s16 m_sScaleFactor;	// Scale factor to which the depth data is downsampled
}T_DataFusion_Set;

typedef struct
{
    f32 b1;
    f32 b2;
    f32 a0;
    f32 a1;
    f32 a2;
    f32 a3;
    f32 coefp;
    f32 coefn;
}T_Gaussian_Rec_St;

typedef struct
{
    T_DataFusion_Set m_Set;
    T_Gaussian_Rec_St m_StGaussian;
    short m_sInitOk;
    s16 m_h;
    s16 m_w;
    s16 m_h_original;
    s16 m_w_original;
    s16 m_nr_shift;
    s16 m_radius;
    u08 m_ucaTable[MLF_APP_2D_RANGE];          // Gaussian on intensity data
    u08 m_ucaTableBF[MLF_APP_INTENSITY_RANGE]; // Gaussian on depth data
    s32 *m_pdaImgBox;
    u16 *m_pt_J_BF_DS;
    s16 m_saDxx[16][16][4];
    //RGB-D Filter
    u08 *m_pt_I_RGB;
    u08 *m_pt_I_GRAY;
    u08 *m_pt_I_R;
    u08 *m_pt_I_G;
    u08 *m_pt_I_B;
    u08 *m_pt_I_GRAY_DS;
    u08 *m_pt_I_R_DS;
    u08 *m_pt_I_G_DS;
    u08 *m_pt_I_B_DS;
    u16 *m_pt_D;
    u16 *m_pt_D_DS;
    u08 *m_pt_Q;
    u08 *m_pt_Q_DS;
    u08 *m_pt_BetaValue;    // Image that contains the beta value to merge J_2 with D
    u08 *m_pt_BetaChannel;  // Image that indicates from which channel (RGB) the Beta Value has been computed
    u16 *m_pt_J;    
    s32 *m_pt_Jk_level0_R;
    s32 *m_pt_Jk_level1_R;
    s32 *m_pt_Jk_level0_G;
    s32 *m_pt_Jk_level1_G;
    s32 *m_pt_Jk_level0_B;
    s32 *m_pt_Jk_level1_B;
    s32 *m_pt_Jk_level0_GRAY;
    s32 *m_pt_Jk_level1_GRAY;
    s32 *m_pt_Wk_R;
    s32 *m_pt_Wk_G;
    s32 *m_pt_Wk_B;
    s32 *m_pt_Wk_GRAY;
#ifdef	MLF_BackGround_Counter
    s32 *m_pt_Jk_level0_R_BG;
    s32 *m_pt_Jk_level1_R_BG;
    s32 *m_pt_Jk_level0_G_BG;
    s32 *m_pt_Jk_level1_G_BG;
    s32 *m_pt_Jk_level0_B_BG;
    s32 *m_pt_Jk_level1_B_BG;
    s32 *m_pt_Jk_level0_GRAY_BG;
    s32 *m_pt_Jk_level1_GRAY_BG;
#endif
    //RGB-D Filter end
}T_DataFusion_St;

/** \brief Data Fusion class.
 * Class that implements low-level data fusion filters to enhance depth maps.
 * \author Frederic GARCIA BECERRO (frederic.garcia@uni.lu)
 */
class c_DataFusion
{
public:
    c_DataFusion(s16 p_sImageWidth=320, s16 p_sImageHeight=240, s16 p_sNr_shift=D_Fusion_ScaleFactor);
    ~c_DataFusion();

    u16* GetEnhancedDepthData(); // Returns the enhanced depth data (after filtering)
    u08* GetGuidanceImageRGB(); // Returns the guidance image I (RGB)
    u08* GetGuidanceImage(); // Returns the guidance image I
    u08* GetGuidanceImage_ds(); // Returns the downsampled guidance image I_ds
    u16* GetDepthMap(); // Returns the input depth map D
    u16* GetDepthMap_ds(); // Returns the downsampled input depth map D_ds
    u08* GetCredibilityMap(); // Returns the credibility map Q
    u08* GetCredibilityMap_ds(); // Returns the downsampled credibility map Q
    u08* GetBlendingMask(); // Returns the blending mask B

    void SetSigmaS(float p_fVal); // Set sigma spatial
    void SetSigmaR(float p_fVal); // Set sigma range
    void SetSigmaQ(float p_fVal); // Set sigma Q

    void DataProcessing(u16* p_usDepthData, u08* p_ucRGBData, s16 p_usFIlterType); // Data processing and internal strcuture update
    /** \brief A Bilateral Filter filter implementation.
     * \note For more information please see
     * <b>C. Tomasi and R. Manduchi. Bilateral filtering for gray and color images.
     * In ICCV, pages 839–846, 1998.</b>
     */
    short BF_Filter(); // Apply Bilateral Filter (BF Filter)
    /** \brief A Joint Bilateral Upsampling (JBU) filter implementation.
     * \note For more information please see
     * <b>J. Kopf, M. Cohen, D. Lischinski, and M. Uyttendaele. Joint Bilateral Upsampling.
     * In SIGGRAPH ’07: ACM SIGGRAPH 2007 papers, page 96, New York, NY, USA, 2007. ACM.</b>
     */
    short JBU_Filter(); // Apply Joint Bilateral Upsampling Filter (JBU Filter)
    /** \brief A Pixel Weighted Average Strategy (PWAS) filter implementation.
     * \note For more information please see
     * <b>F. Garcia, B. Mirbach, B. Ottersten, F. Grandidier, and A. Cuesta. Pixel Weighted Average Strategy for Depth Sensor Data Fusion.
     * In International Conference on Image Processing (ICIP), pages 2805–2808, September 2010.</b>
     */
    short PWAS_Filter();
    /** \brief A Pixel Weighted Average Strategy (PWAS) filter implementation that considers RGB guidance images.
     * \note For more information please see
     * <b>F. Garcia, B. Mirbach, B. Ottersten, F. Grandidier, and A. Cuesta. Pixel Weighted Average Strategy for Depth Sensor Data Fusion.
     * In International Conference on Image Processing (ICIP), pages 2805–2808, September 2010.</b>
     */
    short PWAS_RGB_Filter();
    /** \brief A Unified Multi-Lateral (UML) filter implementation.
     * \note For more information please see
     * <b>F. Garcia, D. Aouada, B. Mirbach, T. Solignac, and B. Ottersten. A New Multi-lateral Filter for Real-Time Depth Enhancement.
     * In Advanced Video and Signal-Based Surveillance (AVSS), 2011.</b>
     */
    short UML_Filter();
    /** \brief A Unified Multi-Lateral (UML) filter implementation  that considers RGB guidance images.
     * \note For more information please see
     * <b>F. Garcia, D. Aouada, B. Mirbach, T. Solignac, and B. Ottersten. A New Multi-lateral Filter for Real-Time Depth Enhancement.
     * In Advanced Video and Signal-Based Surveillance (AVSS), 2011.</b>
     */
    short UML_RGB_Filter();
    /** \brief A Joint Bilateral Upsampling (JBU) filter implementation.
     * \note For more information please see
     * <b>J. Kopf, M. Cohen, D. Lischinski, and M. Uyttendaele. Joint Bilateral Upsampling.
     * In SIGGRAPH ’07: ACM SIGGRAPH 2007 papers, page 96, New York, NY, USA, 2007. ACM.</b>
     */
    short JBU_Filter_Kopf();
    /** \brief A New Joint Bilateral Upsampling (NJBU) filter implementation.
     * \note For more information please see
     * <b>S.-Y. Kim, J.-H. Cho, A. Koschan, and M. Abidi. Spatial and Temporal Enhancement of Depth Images Captured by a Time-of-Flight Depth Sensor.
     * In International Conference on Pattern Recognition (ICPR), pages 2358–2361, August 2010.</b>
     */
    short NJBU_Filter_Kim();
    /** \brief A Noise-Aware Filter for real-time Depth Upsampling filter implementation.
     * \note For more information please see
     * <b>D. Chan, H. Buisman, C. Theobalt, and S. Thrun. A noise-aware filter for real-time depth upsampling.
     * In Workshop on Multi-camera and Multi-modal Sensor Fusion Algorithms and Applications (ECCVW), 2008.</b>
     */
    short NAFDU_Filter_Chan();

private:
    T_DataFusion_St m_DataFusion;
    T_Memory_St m_StMem;
    char m_cMem[D_Fusion_2D_3D_MEM_SIZE]; // Memory handling

    s32 Alloc_Mem(T_Memory_St* p_pSt,char **p_ppVoid, s32 p_lMemSize);
    void Compute_Q_D();
    void Compute_Q_D_RGBD_Filter();
    void ComputeBlendingMask(); // Compute the blending mask (beta) that combines J_2 with J_3 or D
    void ConvertToGrayscale(const unsigned char* p_ptImgIn, long p_lImgSize, unsigned char* p_ptImgOutC);
    void GetMinMaxVal(u08* p_pucMin, u08* p_pucMax, const u08 *p_pucImg, s32 p_lLen);    
    void GetMinMaxVal(u16* p_pusMin, u16* p_pusMax, const u16 *p_pusImg, s32 p_lLen);
    void GetMinMaxVal_RGB(u08* p_pucMin, u08* p_pucMax, const u08 *p_pucImg_R, const u08 *p_pucImg_G, const u08 *p_pucImg_B, s32 p_lLen);
    void Downsample(const u08 *p_pcImgIn, u08 *p_pcImgOut, s16 h, s16 w, s16 scale_exp);
    void Downsample(const u16 *p_pcImgIn, u16 *p_pcImgOut, s16 h, s16 w, s16 scale_exp);
    void Upsample(const u16 *p_psImgIn, u16 *p_psImgOut, s16 h, s16 w, s16 scale_exp);
    void FillColorWeightedTable(u08* p_pucTable, f32 p_pfSigma_range, s16 len);
    void Gaussian_Recursive_Order0_Init(f32 sigma, s16 h, s16 w, T_Gaussian_Rec_St* p_StOut);
    void Gaussian_Recursive_Order0(const T_Gaussian_Rec_St* p_StSet,s32 *p_pdaImg,s32 *p_pdaImgTemp,s16 h,s16 w);
    void Gaussian_Recursive_x(s32 *p_pdImgOut,const s32 *p_pdImgIn, s16 w, s16 h,
                const f32 a0, const f32 p_dA1, const f32 a2, const f32 a3,
                const f32 b1, const f32 b2, const f32 coefp, const f32 coefn);
    void Gaussian_Recursive_y(s32 *p_pdImgOut,const s32 *p_pdImgIn, s16 w, s16 h,
                const f32 a0, const f32 p_dA1, const f32 a2, const f32 a3,
                const f32 b1, const f32 b2, const f32 coefp, const f32 coefn);
};

#endif // C_DATAFUSION_H
