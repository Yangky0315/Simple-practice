#pragma once
#include "RSDIPLib.h"
//#include "d:\study\Ò£¸ÐÍ¼Ïñ´¦Àí\course-rsdip\rsdiplib\include\image.h"
class CImageProcessingEx :
	public CImageProcessing
{
public:
	CImageProcessingEx(void);
	~CImageProcessingEx(void);
	static BOOL histeq(CImageDataset &imgIn, CImageDataset &imgOut);
	static BOOL match(CImageDataset &imgIn1, CImageDataset &imgIn2, CImageDataset &imgOut);
	static BOOL median(int muban,CImageDataset &imgIn,CImageDataset &imgOut);
	static BOOL bilateral(int muban,double s,double r,CImageDataset &imgIn, CImageDataset &imgOut);
	static BOOL laplacian(CImageDataset &imgIn, CImageDataset &imgOut);
	
	static BOOL Creat(CImageDataset &imgOut);
	//static BOOL DFT(CImageDataset &imgIn, CImageDataset &imgOut);
	static BOOL DFT(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, CImageDataset &error);
	static BOOL LG(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, float LD0);
	static BOOL BTWG(CImageDataset &imgIn, CImageDataset &imgOut,CImageDataset &imgOut1, float BD0, int BN);
	static BOOL GG(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, float GD0) ;
	static BOOL QTD(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1);

	static BOOL PCA(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgIn1, CImageDataset &error, int p, int t);
	static BOOL B_I(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, CImageDataset &error);
	//static BOOL I_B(CImageDataset &imgIn, CImageDataset &imgOut);
};

