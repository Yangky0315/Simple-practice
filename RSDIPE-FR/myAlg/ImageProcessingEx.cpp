#include "StdAfx.h"
#include "ImageProcessingEx.h"
#include "math.h"
#include <complex>
#include "gdal_priv.h"
#include "Matrix.h"
const float pi = 3.1415926;
#define MATH_PI 3.141592654


CImageProcessingEx::CImageProcessingEx(void)
{
}

CImageProcessingEx::~CImageProcessingEx(void)
{
}

BOOL CImageProcessingEx::histeq(CImageDataset &imgIn, CImageDataset &imgOut)  
{  
  const int LEVEL = 256;  
  int k, row, col;  
  double hist[LEVEL], sk[LEVEL];  
  /* checking image validation */  
  if(imgIn.empty())  
  {  
    return FALSE;  
  }  
  /* creating the output image */  
  if(FALSE == imgIn.duplicate(imgOut))//get new image from input image  
  {  
    return FALSE;  
  }  
  double *data = imgOut.m_data;  
  /*  步骤 1：计算输入图像的直方图  */  
  for(k=0;k< LEVEL;k++)  //初始化直方图数据为 0  
  {  
    hist[k]=0;
	sk[k]=0;
  }  
  for(row=0; row<imgIn.m_ysize; row++)  //对图像逐行逐列进行扫描，计算完成直方图 
  {  
    for(col=0; col<imgIn.m_xsize; col++)  
    {  
      hist[UINT8(data[row*imgIn.m_xsize+col])]++;  
    }  
  }  
  /*  步骤 2：计算累积直方图  */  
  for(k=1;k< LEVEL-1;k++) //对直方图索引进行扫描，计算完成累积直方图 
  { 
	  sk[k]=sk[k-1]+ hist[k]/ (imgIn.m_ysize*imgIn.m_xsize); 
  } 
  sk[255]=1; 
  /*  步骤 3：根据用累积直方图建立灰度变换函数，并逐行逐列进行像素的灰度变换 */ 
  for(row=0; row<imgOut.m_ysize; row++) 
  { 
	  for(col=0; col<imgOut.m_xsize; col++) 
	  { 
		  for(k=0;k<LEVEL;k++) 
		  { 
			  if(data[row*imgOut.m_xsize+col] == k ) 
			  { 
				  data[row*imgOut.m_xsize+col] = int((LEVEL-1)*sk[k]+0.5);//灰度变换 
				  k=LEVEL; 
			  } 
		  } 
	  } 
  } 
  return TRUE;
}

//直方图匹配
BOOL CImageProcessingEx::match(CImageDataset &imgIn1, CImageDataset &imgIn2,CImageDataset &imgOut)  
{
	const int LEVEL = 256; 
	int k, row, col,m,s; 
	double z,min,hist_1[LEVEL],hist_2[LEVEL], c_hist_1[LEVEL],c_hist_2[LEVEL],c_hist,m_1[LEVEL]; 
	/* checking image validation */ 
	if(imgIn1.empty()) 
	{ 
		return FALSE; 
	} 
	/* creating the output image */ 
	if(FALSE == imgIn1.duplicate(imgOut))//get new image from input image 
	{ 
		return FALSE; 
	} 
	double *data = imgOut.m_data; 
	/* 步骤 1：计算输入图像 1 的直方图和累积直方图 */ 
	for(k=0;k< LEVEL;k++) 
		//初始化直方图数据为 0 
	{ 
		hist_1[k]=0; 
		c_hist_1[k]=0;
		c_hist_2[k]=0;
		m_1[k]=0;
	} 
	for(row=0; row<imgIn1.m_ysize; row++) //对图像逐行逐列进行扫描，计算完成直方图 
	{
		for(col=0; col<imgIn1.m_xsize; col++) 
		{ 
			hist_1[UINT8(data[row*imgIn1.m_xsize+col])]++; 
		} 
	} 
	c_hist_1[0]= hist_1[0]/ (imgIn1.m_ysize*imgIn1.m_xsize); 
	for(k=1;k< LEVEL-1;k++) //对直方图索引进行扫描，计算完成累积直方图 
	{ 
		c_hist_1[k]=c_hist_1[k-1]+ hist_1[k]/(imgIn1.m_ysize*imgIn1.m_xsize); 
	} 
	c_hist_1[255]=1; 
	/* 步骤 2：计算输入图像 2 的直方图和累积直方图 */ 
	for(k=0;k< LEVEL;k++) 
		//初始化直方图数据为 0 
	{ 
		hist_2[k]=0; 
	} 
	for(row=0; row<imgIn2.m_ysize; row++) //对图像逐行逐列进行扫描，计算完成直方图 
	{ 
		for(col=0; col<imgIn2.m_xsize; col++) 
		{ 
			hist_2[UINT8(imgIn2.m_data[row*imgIn2.m_xsize+col])]++; 
		} 
	} 
	c_hist_2[0]= hist_2[0]/ (imgIn2.m_ysize*imgIn2.m_xsize); 
	for(k=1;k< LEVEL-1;k++) //对直方图索引进行扫描，计算完成累积直方图 
	{ 
		c_hist_2[k]=c_hist_2[k-1]+ hist_2[k]/ (imgIn2.m_ysize*imgIn2.m_xsize); 
	} 
	c_hist_2[255]=1; 
	/* 步骤 3：进行匹配 */ 
	min=2; 
	z=0; 
	s=0; 
	for(k=0;k<LEVEL;k++) 
	{ 
		if(c_hist_2[k]==0) 
			z=z+1; 
	} 
	for(k=0;k<LEVEL;k++) 
	{ 
		min=2;
		for(m=z;m<LEVEL;m++) 
		{ 
			c_hist=fabs(c_hist_1[k]-c_hist_2[m]); 
			if(c_hist<min) 
			{ 
				min=c_hist; 
				s=m; 
			} 
		} 
		m_1[k]=s; 
	} 
	for(row=0; row<imgOut.m_ysize; row++) 
	{ 
		for(col=0; col<imgOut.m_xsize; col++) 
		{ 
			for(k=0;k<LEVEL;k++) 
			{ 
				if(data[row*imgOut.m_xsize+col] == k ) 
				{ 
					data[row*imgOut.m_xsize+col] = m_1[k];//灰度变换 
					k=LEVEL; 
				} 
			} 
		} 
	} 
	return TRUE; 

}

//中值滤波
BOOL CImageProcessingEx::median(int muban,CImageDataset &imgIn,CImageDataset &imgOut) 
{ 
	int k, row, col,i,j; 
	double z,s[10000]; 
	/* checking image validation */ 
	if(imgIn.empty()) 
	{ 
		return FALSE; 
	} 
	/* creating the output image */ 
	if(FALSE == imgIn.duplicate(imgOut))//get new image from input image 
	{ 
		return FALSE; 
	} 
	double *data = imgOut.m_data; 
	k=0; 
	for(k=0;k<1000;k++) //初始化数据为 0 
	{ 
		s[k]=0; 
	} 
	k=0;
	for(row=(0+(muban-1)/2); row<(imgOut.m_ysize-(muban-1)/2); row++)
	{
		for(col=(0+(muban-1)/2); col<(imgOut.m_xsize-(muban-1)/2); col++) 
		{ 
			k=0; 
			for(i=0;i<muban;i++) 
			{ 
				for(j=0;j<muban;j++) 
				{ 
					s[k]=data[(row-(muban- 1)/2+i)*(imgOut.m_xsize)+col+j-(muban-1)/2] ;//灰度变换 
					k=k+1; 
				} 
			} 
			for(i=0;i<k-1;i++)//冒泡排序计算 
			{ 
				for(j=i;j<(k-1-i);j++) 
				{ 
					if(s[j]<s[j+1]) 
					{ 
						z=s[j]; 
						s[j]=s[j+1]; 
						s[j+1]=z; 
					} 
				} 
			} 
			data[row*imgOut.m_xsize+col]=s[(k-1)/2];//计算新的灰度值 
		} 
	} 
	return TRUE;
}

//双边滤波
BOOL CImageProcessingEx::bilateral(int muban,double s,double r,CImageDataset &imgIn, CImageDataset &imgOut) 
{ 
	int k, row, col,i,j; 
	double s_1,s_2,s_3,s_4; 
	/* checking image validation */ 
	if(imgIn.empty()) 
	{ 
		return FALSE; 
	} 
	/* creating the output image */ 
	if(FALSE == imgIn.duplicate(imgOut))//get new image from input image 
	{ 
		return FALSE; 
	} 
	double *data = imgIn.m_data; 
	for(row=0; row<imgOut.m_ysize; row++) 
	{ 
		for(col=0; col<imgOut.m_xsize; col++) 
		{ 
			data[row*imgOut.m_xsize+col]=data[row*imgOut.m_xsize+col]/255;//归一化灰度 
		} 
	} 
	for(row=(0+(muban-1)/2); row<(imgOut.m_ysize-(muban-1)/2); row++) 
	{ 
		for(col=(0+(muban-1)/2); col<(imgOut.m_xsize-(muban-1)/2); col++) 
		{ 
			s_3=0.0; 
			s_4=0.0; 
			k=0; 
			for(i=0;i<muban;i++) 
			{ 
				for(j=0;j<muban;j++) 
				{ 
					s_1=exp((-1.0*(i+1)*(i+1)+(j*j)*(j+1))/(2*s*s)); 
					s_2=exp(-pow((data[row*imgOut.m_xsize+col]- data[(row-(muban-1)/2+i)*(imgOut.m_xsize)+col+j-(muban-1)/2]),2)/(2*r*r));
					s_3=s_3+s_1*s_2; 
					s_4=s_4+s_1*s_2*data[(row-(muban- 
						1)/2+i)*(imgOut.m_xsize)+col+j-(muban-1)/2]; 
				} 
			} 
			imgOut.m_data[row*imgOut.m_xsize+col]=(s_4/s_3)*255;//灰度变换 
		} 
	} 
	return TRUE;
}

//锐化处理
BOOL CImageProcessingEx::laplacian(CImageDataset &imgIn, CImageDataset &imgOut)  
{
	int k, row, col,i,j,m;
	double s_1;
	/* checking image validation */
	if(imgIn.empty())
	{
		return FALSE;
	}
	/* creating the output image */
	if(FALSE == imgIn.duplicate(imgOut))//get new image from input image
	{
		return FALSE;
	}
	double *data = imgIn.m_data;
	m=3;
	int m_x[9]={0,1,0,1,-4,1,0,1,0};//定义拉普拉斯模板
	for(row=(0+(m-1)/2); row<(imgOut.m_ysize-(m-1)/2); row++)
	{
		for(col=(0+(m-1)/2); col<(imgOut.m_xsize-(m-1)/2); col++)
		{
			s_1=0.0;
			k=0;
			for(i=0;i<m;i++)
			{
				for(j=0;j<m;j++)
				{
					s_1=s_1+data[(row-(m-1)/2+i)*(imgOut.m_xsize)+col+j-(m-1)/2]*m_x[k];//卷积运算
					k=k+1;
				} 
			}
			imgOut.m_data[row*imgOut.m_xsize+col]=data[row*imgOut.m_xsize+col]-s_1;//得到锐化结果
		} 
	}
	return TRUE;
}

//创建测试图像
BOOL CImageProcessingEx::Creat(CImageDataset &imgOut) 
{
	imgOut.m_xsize = 512;
	imgOut.m_ysize = 512;
	imgOut.m_rastercount = 1;
	imgOut.m_data = (double*)CPLMalloc(sizeof(double)*imgOut.m_ysize*imgOut.m_xsize*imgOut.m_rastercount);
	//图像宽20，高40
	int x, y;
	for (y=0; y<imgOut.m_ysize; y++)                  
	{
		for (x=0; x<imgOut.m_xsize; x++)
		{
			if ((x<=266 && x>=246) && (y<=276 && y>=236))
				imgOut.m_data[y*imgOut.m_xsize+x]=255;
			else
				imgOut.m_data[y*imgOut.m_xsize+x]=0;
		}
	}

	return TRUE;
}

//二维离散傅里叶
//BOOL CImageProcessingEx::DFT(CImageDataset &imgIn, CImageDataset &imgOut)
//{
//	///* checking image validation */
//	//if(imgIn.empty())
//	//{
//	//	return FALSE;
//	//}
//
//	///* prepare output image */
//	//if(FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 5))// fourier, angle, energy, real, imagery
//	//{
//	//	return FALSE;
//	//}
//
//	///* prepare temporary image */
//	//CImageDataset imgTmp;
//	//if(FALSE == imgIn.duplicate(imgTmp))//get new image from input image
//	//{
//	//	return FALSE;
//	//}
//	//for(int row=0; row<imgTmp.m_ysize; row++)
//	//{
//	//	for(int col=0; col<imgTmp.m_xsize; col++)
//	//	{
//	//		if ((row+col)%2==1) 
//	//		{
//	//			imgTmp.m_data[row*imgIn.m_xsize+col]*=-1 ;
//	//		}
//	//	}
//	//}
//
//	//using namespace std;
//	//complex <double> *fxv=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
//	//complex <double> *fuv=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
//	//for(int x=0;x<imgTmp.m_ysize;x++)
//	//{
//	//	for(int v=0;v<imgTmp.m_xsize;v++)
//	//	{
//	//		fxv[x*imgTmp.m_xsize+v] = complex <double> (0,0);
//	//		for(int y=0; y<imgTmp.m_xsize;y++)
//	//		{
//	//			complex <double> comp(0,-2*MATH_PI*v*y/imgTmp.m_xsize);
//	//			fxv[x*imgTmp.m_xsize+v] += imgTmp.m_data[x*imgTmp.m_xsize+y]*exp(comp);
//	//		}
//	//		fxv[x*imgTmp.m_xsize+v] /= imgTmp.m_xsize;
//	//	}
//	//}
//	//for(int v=0; v<imgTmp.m_xsize;v++)
//	//{
//	//	for(int u=0; u<imgTmp.m_ysize; u++)
//	//	{
//	//		fuv[u*imgTmp.m_xsize+v] = complex <double> (0,0);
//	//		for(int x=0; x<imgTmp.m_ysize; x++)
//	//		{
//	//			complex <double> comp(0,-2*MATH_PI*u*x/imgTmp.m_ysize);
//	//			fuv[u*imgTmp.m_xsize+v] += fxv[x*imgTmp.m_xsize+v]*exp(comp);
//	//		}
//	//		fuv[u*imgTmp.m_xsize+v] /= imgTmp.m_ysize;
//	//	}
//	//}
//
//	//double *data = imgOut.m_data;
//	//for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
//	//{
//	//	data[p] = log(norm(fuv[p])+1);//傅立叶谱
//	//	data[p+imgTmp.m_xsize*imgTmp.m_ysize] = arg(fuv[p]);//角谱
//	//	data[p+imgTmp.m_xsize*imgTmp.m_ysize*2]=norm(fuv[p]);//能量谱
//	//	data[p+imgTmp.m_xsize*imgTmp.m_ysize*3]=fuv[p].real();//实部
//	//	data[p+imgTmp.m_xsize*imgTmp.m_ysize*4]=fuv[p].imag();//虚部
//	//}
//
//	//delete fxv;
//	//delete fuv;
//
//	return TRUE;
//}


//BOOL CImageProcessingEx::IDFT(CImageDataset &imgIn, CImageDataset &imgOut)
//{
//	/* checking image validation */
//	if(imgIn.empty())
//	{
//		return FALSE;
//	}
//	CImageDataset imgTmp;
//	if(FALSE == imgIn.duplicate(imgTmp))//get new image from input image
//	{
//		return FALSE;
//	}
//	int row, col;
//	for(row=0; row<imgTmp.m_ysize; row++)
//	{
//		for(col=0; col<imgTmp.m_xsize; col++)
//		{
//			if ((row+col)%2==1) 
//			{
//				imgTmp.m_data[row*imgIn.m_xsize+col]*=-1 ;
//			} } }
//	using namespace std;
//	complex <double> *fxv=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
//	complex <double> *fuv=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
//	complex <double> *fuy=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
//	complex <double> *fxy=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
//	//对输入图像每一行做一维 DFT
//	for(int x=0;x<imgTmp.m_ysize;x++)
//	{
//		for(int v=0;v<imgTmp.m_xsize;v++)
//		{
//			fxv[x*imgTmp.m_xsize+v] = complex <double> (0,0);
//			for(int y=0; y<imgTmp.m_xsize;y++)
//			{
//				complex <double> comp(0,-2*MATH_PI*v*y/imgTmp.m_xsize);
//				fxv[x*imgTmp.m_xsize+v] += 
//					imgTmp.m_data[x*imgTmp.m_xsize+y]*exp(comp);
//			}
//			fxv[x*imgTmp.m_xsize+v] /= imgTmp.m_xsize;
//		} }
//	//对矩阵 Fxv 的每一列做 DFT
//	for(int v=0; v<imgTmp.m_xsize;v++)
//	{
//		for(int u=0; u<imgTmp.m_ysize; u++)
//		{
//			fuv[u*imgTmp.m_xsize+v] = complex <double> (0,0);
//			for(int x=0; x<imgTmp.m_ysize; x++)
//			{
//				complex <double> comp(0,-2*MATH_PI*u*x/imgTmp.m_ysize);
//				fuv[u*imgTmp.m_xsize+v] += fxv[x*imgTmp.m_xsize+v]*exp(comp);
//			}
//			fuv[u*imgTmp.m_xsize+v] /= imgTmp.m_ysize;
//		} }
//	//对图像每一行做一维 IDFT
//	for(int u=0;u<imgTmp.m_ysize;u++)
//	{
//		for(int y=0;y<imgTmp.m_xsize;y++)
//		{
//			fuy[u*imgTmp.m_xsize+y] =complex <double> (0,0);
//			for(int y1=0; y1<imgTmp.m_xsize;y1++)
//			{
//				complex <double> comp(0,2*MATH_PI*y*y1/imgTmp.m_xsize);
//				fuy[u*imgTmp.m_xsize+y] += fuv[u*imgTmp.m_xsize+y1]*exp(comp);
//			}
//			//fvy[v*imgTmp.m_xsize+y] /= imgTmp.m_xsize;
//		} }
//	//对矩阵 Fvy 的每一列做 IDFT
//	for(int y=0; y<imgTmp.m_xsize;y++)
//	{
//		for(int x=0; x<imgTmp.m_ysize; x++)
//		{
//			fxy[x*imgTmp.m_xsize+y] =complex <double> (0,0);
//			for(int x1=0; x1<imgTmp.m_ysize; x1++)
//			{
//				complex <double> comp(0,2*MATH_PI*x*x1/imgTmp.m_ysize);
//				fxy[x*imgTmp.m_xsize+y] += fuy[x1*imgTmp.m_xsize+y]*exp(comp);
//			}
//			//fxy[x*imgTmp.m_xsize+y] /= imgTmp.m_ysize;
//		} }if(FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 2))// DIFT,ERROR
//		{
//			return FALSE;
//		}
//		double *data = imgOut.m_data;
//		for(row=0; row<imgOut.m_ysize; row++)
//		{
//			for(col=0; col<imgOut.m_xsize; col++)
//			{
//				if ((row+col)%2==1) 
//				{
//					data[row*imgOut.m_xsize+col] = -1*fxy[row*imgOut.m_xsize+col].real(); 
//					//DIFT
//				}
//				else
//				{
//					data[row*imgOut.m_xsize+col] = fxy[row*imgOut.m_xsize+col].real();
//				} } }
//		//计算离散傅里叶变换的结果与测试图像之间的误差
//		for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
//		{
//			data[p+imgTmp.m_xsize*imgTmp.m_ysize]=abs(imgIn.m_data[p]-data[p]);
//		}
//		delete fxv;
//		delete fuv;
//		delete fuy;
//		delete fxy;
//		return true;
//
//}


BOOL CImageProcessingEx::DFT(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, CImageDataset &error)
{
	/* checking image validation */
	if(imgIn.empty())
	{
		return FALSE;
	}
		/* prepare output image */
	if(FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 5))// fourier, angle, energy, real, imagery
	{
		return FALSE;
	}
		/* prepare output image */
	if(FALSE == imgOut1.create(imgIn.m_xsize, imgIn.m_ysize, 5))// fourier, angle, energy, real, imagery
	{
		return FALSE;
	}
	if(FALSE == error.create(imgIn.m_xsize, imgIn.m_ysize, 5))// fourier, angle, energy, real, imagery
	{
		return FALSE;
	}
	CImageDataset imgTmp;
	if(FALSE == imgIn.duplicate(imgTmp))//get new image from input image
	{
		return FALSE;
	}
	int row, col;
	for(row=0; row<imgTmp.m_ysize; row++)
	{
		for(col=0; col<imgTmp.m_xsize; col++)
		{
			if ((row+col)%2==1) 
			{
				imgTmp.m_data[row*imgIn.m_xsize+col]*=-1 ;
			} 
		} 
	}
	using namespace std;
	complex <double> *fxv=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuv=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuy=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fxy=new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	//对输入图像每一行做一维 DFT
	for(int x=0;x<imgTmp.m_ysize;x++)
	{
		for(int v=0;v<imgTmp.m_xsize;v++)
		{
			fxv[x*imgTmp.m_xsize+v] = complex <double> (0,0);
			for(int y=0; y<imgTmp.m_xsize;y++)
			{
				complex <double> comp(0,-2*MATH_PI*v*y/imgTmp.m_xsize);
				fxv[x*imgTmp.m_xsize+v] += 
					imgTmp.m_data[x*imgTmp.m_xsize+y]*exp(comp);
			}
			fxv[x*imgTmp.m_xsize+v] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fxv 的每一列做 DFT
	for(int v=0; v<imgTmp.m_xsize;v++)
	{
		for(int u=0; u<imgTmp.m_ysize; u++)
		{
			fuv[u*imgTmp.m_xsize+v] = complex <double> (0,0);
			for(int x=0; x<imgTmp.m_ysize; x++)
			{
				complex <double> comp(0,-2*MATH_PI*u*x/imgTmp.m_ysize);
				fuv[u*imgTmp.m_xsize+v] += fxv[x*imgTmp.m_xsize+v]*exp(comp);
			}
			fuv[u*imgTmp.m_xsize+v] /= imgTmp.m_ysize;
		} 
	}

    double *data_1 = imgOut1.m_data;
	for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
	{
		data_1[p] = abs(fuv[p]);//取绝对值输出
	}

	//对图像每一行做一维 IDFT
	for(int u=0;u<imgTmp.m_ysize;u++)
	{
		for(int y=0;y<imgTmp.m_xsize;y++)
		{
			fuy[u*imgTmp.m_xsize+y] =complex <double> (0,0);
			for(int y1=0; y1<imgTmp.m_xsize;y1++)
			{
				complex <double> comp(0,2*MATH_PI*y*y1/imgTmp.m_xsize);
				fuy[u*imgTmp.m_xsize+y] += fuv[u*imgTmp.m_xsize+y1]*exp(comp);
			}
			//fvy[v*imgTmp.m_xsize+y] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fvy 的每一列做 IDFT
	for(int y=0; y<imgTmp.m_xsize;y++)
	{
		for(int x=0; x<imgTmp.m_ysize; x++)
		{
			fxy[x*imgTmp.m_xsize+y] =complex <double> (0,0);
			for(int x1=0; x1<imgTmp.m_ysize; x1++)
			{
				complex <double> comp(0,2*MATH_PI*x*x1/imgTmp.m_ysize);
				fxy[x*imgTmp.m_xsize+y] += fuy[x1*imgTmp.m_xsize+y]*exp(comp);
			}
			//fxy[x*imgTmp.m_xsize+y] /= imgTmp.m_ysize;
		} 
	}
	//if(FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 2))// DIFT,ERROR
	//{
	//	return FALSE;
	//}
	double *data = imgOut.m_data;
	for(row=0; row<imgOut.m_ysize; row++)
	{
		for(col=0; col<imgOut.m_xsize; col++)
		{
			if ((row+col)%2==1) 
			{
				data[row*imgOut.m_xsize+col] = -1*fxy[row*imgOut.m_xsize+col].real(); 
				//DIFT
			}
			else
			{
				data[row*imgOut.m_xsize+col] = fxy[row*imgOut.m_xsize+col].real();
			} 
		} 
	}
	//计算离散傅里叶变换的结果与测试图像之间的误差	
	double *data_error = error.m_data;
	for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
	{
		data_error[p+imgTmp.m_xsize*imgTmp.m_ysize]=abs(imgIn.m_data[p]-data[p]);
	}
	delete fxv;
	delete fuv;
	delete fuy;
	delete fxy;
	return TRUE;
}

BOOL CImageProcessingEx::LG(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, float D0)
{
	if (imgIn.empty())
	{
		return FALSE;
	}
	CImageDataset imgTmp;
	if (FALSE == imgIn.duplicate(imgTmp))//get new image from input image
	{
		return FALSE;
	}
	int row, col;
	for (row = 0; row < imgTmp.m_ysize; row++)
	{
		for (col = 0; col < imgTmp.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				imgTmp.m_data[row*imgIn.m_xsize + col] *= -1;
			} 
		} 
	}
	using namespace std;
	complex <double> *fxv = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuv = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuy = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fxy = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	//对输入图像每一行做一维 DFT
	for (int x = 0; x < imgTmp.m_ysize; x++)
	{
		for (int v = 0; v < imgTmp.m_xsize; v++)
		{
			fxv[x*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int y = 0; y < imgTmp.m_xsize; y++)
			{
				complex <double> comp(0, -2 * MATH_PI*v*y / imgTmp.m_xsize);
				fxv[x*imgTmp.m_xsize + v] += imgTmp.m_data[x*imgTmp.m_xsize + y] * 
					exp(comp);
			}
			fxv[x*imgTmp.m_xsize + v] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fxv 的每一列做 DFT
	for (int v = 0; v < imgTmp.m_xsize; v++)
	{
		for (int u = 0; u < imgTmp.m_ysize; u++)
		{
			fuv[u*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int x = 0; x < imgTmp.m_ysize; x++)
			{
				complex <double> comp(0, -2 * MATH_PI*u*x / imgTmp.m_ysize);
				fuv[u*imgTmp.m_xsize + v] += fxv[x*imgTmp.m_xsize + v] * exp(comp);
			}
			fuv[u*imgTmp.m_xsize + v] /= imgTmp.m_ysize;
		} 
	}
	if (FALSE == imgOut1.create(imgIn.m_xsize, imgIn.m_ysize, 3))// DIFT,ERROR
	{
		return FALSE;
	}
	double *data_1 = imgOut1.m_data;
	for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
	{
		data_1[p] = abs(fuv[p]);//傅立叶谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize] = arg(fuv[p]);//角谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize*2]=norm(fuv[p]);//能量谱
	}
	//理想高通滤波器
	for (row = 0; row < imgIn.m_ysize; row++)
	{
		for (col = 0; col < imgIn.m_xsize; col++)
		{
			// 阻止低频信息通过
			if (sqrt(pow(col - imgIn.m_xsize / 2.0, 2) + pow(row - imgIn.m_ysize / 2.0, 2)) <= D0)
			{
				fuv[row*imgTmp.m_xsize + col] = complex <double>(0, 0);
			} 
		} 
	}
	//对图像每一行做一维 IDFT
	for (int u = 0; u < imgTmp.m_ysize; u++)
	{
		for (int y = 0; y < imgTmp.m_xsize; y++)
		{
			fuy[u*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int y1 = 0; y1 < imgTmp.m_xsize; y1++)
			{
				complex <double> comp(0, 2 * MATH_PI*y*y1 / imgTmp.m_xsize);
				fuy[u*imgTmp.m_xsize + y] += fuv[u*imgTmp.m_xsize + y1] * exp(comp);
			}
			//fvy[v*imgTmp.m_xsize+y] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fvy 的每一列做 IDFT
	for (int y = 0; y < imgTmp.m_xsize; y++)
	{
		for (int x = 0; x < imgTmp.m_ysize; x++)
		{
			fxy[x*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int x1 = 0; x1 < imgTmp.m_ysize; x1++)
			{
				complex <double> comp(0, 2 * MATH_PI*x*x1 / imgTmp.m_ysize);
				fxy[x*imgTmp.m_xsize + y] += fuy[x1*imgTmp.m_xsize + y] * exp(comp);
			}
			//fxy[x*imgTmp.m_xsize+y] /= imgTmp.m_ysize;
		} 
	}
	if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 2))// DIFT,ERROR
	{
		return FALSE;
	}
	double *data = imgOut.m_data;
	for (row = 0; row < imgOut.m_ysize; row++)
	{
		for (col = 0; col < imgOut.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				data[row*imgOut.m_xsize + col] = -1 * fxy[row*imgOut.m_xsize + 
					col].real(); //DIFT
			}
			else
			{
				data[row*imgOut.m_xsize + col] = fxy[row*imgOut.m_xsize + col].real();
			} 
		} 
	}
	//计算离散傅里叶变换的结果与测试图像之间的误差
	for (int p = 0; p < imgTmp.m_xsize*imgTmp.m_ysize; p++)
	{
		data[p + imgTmp.m_xsize*imgTmp.m_ysize] = abs(imgIn.m_data[p] - data[p]);
	}
	delete fxv;
	delete fuv;
	delete fuy;
	delete fxy;
	return TRUE;
}

BOOL CImageProcessingEx::GG(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, float D0) {
	/* checking image validation */
	if (imgIn.empty())
	{
		return FALSE;
	}
	CImageDataset imgTmp;
	if (FALSE == imgIn.duplicate(imgTmp))//get new image from input image
	{
		return FALSE;
	}
	int row, col;
	for (row = 0; row < imgTmp.m_ysize; row++)
	{
		for (col = 0; col < imgTmp.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				imgTmp.m_data[row*imgIn.m_xsize + col] *= -1;
			} 
		} 
	}
	using namespace std;
	complex <double> *fxv = new complex 
		<double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuv = new complex 
		<double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuy = new complex 
		<double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fxy = new complex 
		<double>[imgIn.m_xsize*imgIn.m_ysize];
	//对输入图像每一行做一维 DFT
	for (int x = 0; x < imgTmp.m_ysize; x++)
	{
		for (int v = 0; v < imgTmp.m_xsize; v++)
		{
			fxv[x*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int y = 0; y < imgTmp.m_xsize; y++)
			{
				complex <double> comp(0, -2 * MATH_PI*v*y / imgTmp.m_xsize);
				fxv[x*imgTmp.m_xsize + v] += 
					imgTmp.m_data[x*imgTmp.m_xsize + y] * exp(comp);
			}
			fxv[x*imgTmp.m_xsize + v] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fxv 的每一列做 DFT
	for (int v = 0; v < imgTmp.m_xsize; v++)
	{
		for (int u = 0; u < imgTmp.m_ysize; u++)
		{
			fuv[u*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int x = 0; x < imgTmp.m_ysize; x++)
			{
				complex <double> comp(0, -2 * MATH_PI*u*x / imgTmp.m_ysize);
				fuv[u*imgTmp.m_xsize + v] += fxv[x*imgTmp.m_xsize + v] * 
					exp(comp);
			}
			fuv[u*imgTmp.m_xsize + v] /= imgTmp.m_ysize;
		} 
	}
	if (FALSE == imgOut1.create(imgIn.m_xsize, imgIn.m_ysize, 3))// DIFT,ERROR
	{
		return FALSE;
	}
	double *data_1 = imgOut1.m_data;
	for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
	{
		data_1[p] = abs(fuv[p]);//傅立叶谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize] = arg(fuv[p]);//角谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize*2]=norm(fuv[p]);//能量谱
	}
	//高斯高通滤波器
	for (row = 0; row < imgIn.m_ysize; row++)
	{
		for (col = 0; col < imgIn.m_xsize; col++)
		{
			// 构建高斯滤波器的传递函数
			double D = sqrt(pow(col - imgIn.m_xsize / 2.0, 2) + pow(row - imgIn.m_ysize / 2.0, 2));
			fuv[row*imgTmp.m_xsize + col] = fuv[row*imgTmp.m_xsize + col] * (1 - exp(-0.5*pow(D / D0, 2)));
		} 
	}
	//对图像每一行做一维 IDFT
	for (int u = 0; u < imgTmp.m_ysize; u++)
	{
		for (int y = 0; y < imgTmp.m_xsize; y++)
		{
			fuy[u*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int y1 = 0; y1 < imgTmp.m_xsize; y1++)
			{
				complex <double> comp(0, 2 * MATH_PI*y*y1 / 
					imgTmp.m_xsize);
				fuy[u*imgTmp.m_xsize + y] += fuv[u*imgTmp.m_xsize + y1] * 
					exp(comp);
			}
		} 
	}
	//对矩阵 Fvy 的每一列做 IDFT
	for (int y = 0; y < imgTmp.m_xsize; y++)
	{
		for (int x = 0; x < imgTmp.m_ysize; x++)
		{
			fxy[x*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int x1 = 0; x1 < imgTmp.m_ysize; x1++)
			{
				complex <double> comp(0, 2 * MATH_PI*x*x1 / 
					imgTmp.m_ysize);
				fxy[x*imgTmp.m_xsize + y] += fuy[x1*imgTmp.m_xsize + y] * exp(comp);
			}
			//fxy[x*imgTmp.m_xsize+y] /= imgTmp.m_ysize;
		} 
	}
	if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 2))// DIFT,ERROR
	{
		return FALSE;
	}
	double *data = imgOut.m_data;
	for (row = 0; row < imgOut.m_ysize; row++)
	{
		for (col = 0; col < imgOut.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				data[row*imgOut.m_xsize + col] = -1 * fxy[row*imgOut.m_xsize + col].real(); //DIFT
			}
			else
			{
				data[row*imgOut.m_xsize + col] = fxy[row*imgOut.m_xsize + col].real();
			} 
		} 
	}
	//计算离散傅里叶变换的结果与测试图像之间的误差
	for (int p = 0; p < imgTmp.m_xsize*imgTmp.m_ysize; p++)
	{
		data[p + imgTmp.m_xsize*imgTmp.m_ysize] = abs(imgIn.m_data[p] -	data[p]);
	}
	delete fxv;
	delete fuv;
	delete fuy;
	delete fxy;
	return TRUE;
}

BOOL CImageProcessingEx::BTWG(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, float D0,int n) {
	/* checking image validation */
	if (imgIn.empty())
	{
		return FALSE;
	}
	CImageDataset imgTmp;
	if (FALSE == imgIn.duplicate(imgTmp))//get new image from input image
	{
		return FALSE;
	}
	int row, col;
	for (row = 0; row < imgTmp.m_ysize; row++)
	{
		for (col = 0; col < imgTmp.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				imgTmp.m_data[row*imgIn.m_xsize + col] *= -1;
			} 
		}
	}
	using namespace std;
	complex <double> *fxv = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuv = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuy = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fxy = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	//对输入图像每一行做一维 DFT
	for (int x = 0; x < imgTmp.m_ysize; x++)
	{
		for (int v = 0; v < imgTmp.m_xsize; v++)
		{
			fxv[x*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int y = 0; y < imgTmp.m_xsize; y++) {
				complex <double> comp(0, -2 * MATH_PI*v*y / imgTmp.m_xsize);
				fxv[x*imgTmp.m_xsize + v] += imgTmp.m_data[x*imgTmp.m_xsize + y] * 
					exp(comp);
			}
			fxv[x*imgTmp.m_xsize + v] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fxv 的每一列做 DFT
	for (int v = 0; v < imgTmp.m_xsize; v++)
	{
		for (int u = 0; u < imgTmp.m_ysize; u++)
		{
			fuv[u*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int x = 0; x < imgTmp.m_ysize; x++)
			{
				complex <double> comp(0, -2 * MATH_PI*u*x / imgTmp.m_ysize);
				fuv[u*imgTmp.m_xsize + v] += fxv[x*imgTmp.m_xsize + v] * exp(comp);
			}
			fuv[u*imgTmp.m_xsize + v] /= imgTmp.m_ysize;
		} 
	}
	if (FALSE == imgOut1.create(imgIn.m_xsize, imgIn.m_ysize, 3))// DIFT,ERROR
	{
		return FALSE;
	}
	double *data_1 = imgOut1.m_data;
	for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
	{
		data_1[p] = abs(fuv[p]);//傅立叶谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize] = arg(fuv[p]);//角谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize*2]=norm(fuv[p]);//能量谱
	}
	//ButterWorth 滤波器
	for (row = 0; row < imgIn.m_ysize; row++)
	{
		for (col = 0; col < imgIn.m_xsize; col++)
		{
			//构建 n 阶 ButterWorth 高通滤波器传递函数
			double D = sqrt(pow(col - imgIn.m_xsize / 2.0, 2) + pow(row - imgIn.m_ysize / 2.0, 2));
			fuv[row*imgTmp.m_xsize + col] = fuv[row*imgTmp.m_xsize + col] * (1 / (1 + pow(D0 / D, n)));
		} 
	}
	//对图像每一行做一维 IDFT
	for (int u = 0; u < imgTmp.m_ysize; u++)
	{
		for (int y = 0; y < imgTmp.m_xsize; y++)
		{
			fuy[u*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int y1 = 0; y1 < imgTmp.m_xsize; y1++)
			{
				complex <double> comp(0, 2 * MATH_PI*y*y1 / imgTmp.m_xsize);
				fuy[u*imgTmp.m_xsize + y] += fuv[u*imgTmp.m_xsize + y1] * exp(comp);
			}
			//fvy[v*imgTmp.m_xsize+y] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fvy 的每一列做 IDFT
	for (int y = 0; y < imgTmp.m_xsize; y++)
	{
		for (int x = 0; x < imgTmp.m_ysize; x++)
		{
			fxy[x*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int x1 = 0; x1 < imgTmp.m_ysize; x1++)
			{
				complex <double> comp(0, 2 * MATH_PI*x*x1 / imgTmp.m_ysize);
				fxy[x*imgTmp.m_xsize + y] += fuy[x1*imgTmp.m_xsize + y] * exp(comp);
			}
			//fxy[x*imgTmp.m_xsize+y] /= imgTmp.m_ysize;
		} 
	}
	if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 2))// DIFT,ERROR
	{
		return FALSE;
	}
	double *data = imgOut.m_data;
	for (row = 0; row < imgOut.m_ysize; row++)
	{
		for (col = 0; col < imgOut.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				data[row*imgOut.m_xsize + col] = -1 * fxy[row*imgOut.m_xsize + 
					col].real(); //DIFT
			}
			else
			{
				data[row*imgOut.m_xsize + col] = fxy[row*imgOut.m_xsize + col].real();
			} 
		} 
	}
	//计算离散傅里叶变换的结果与测试图像之间的误差
	for (int p = 0; p < imgTmp.m_xsize*imgTmp.m_ysize; p++)
	{
		data[p + imgTmp.m_xsize*imgTmp.m_ysize] = abs(imgIn.m_data[p] - data[p]);
	}
	delete fxv;
	delete fuv;
	delete fuy;
	delete fxy;
	return TRUE;
}

BOOL CImageProcessingEx::QTD(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1) {
	/* checking image validation */
	if (imgIn.empty())
	{
		return FALSE;
	}
	CImageDataset imgTmp;
	if (FALSE == imgIn.duplicate(imgTmp))//get new image from input image
	{
		return FALSE;
	}
	int row, col;
	for (row = 0; row < imgTmp.m_ysize; row++)
	{
		for (col = 0; col < imgTmp.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				imgTmp.m_data[row*imgIn.m_xsize + col] *= -1;
			} 
		} 
	}
	using namespace std;
	complex <double> *fxv = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuv = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fuy = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	complex <double> *fxy = new complex <double>[imgIn.m_xsize*imgIn.m_ysize];
	//对输入图像每一行做一维 DFT
	for (int x = 0; x < imgTmp.m_ysize; x++)
	{
		for (int v = 0; v < imgTmp.m_xsize; v++)
		{
			fxv[x*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int y = 0; y < imgTmp.m_xsize; y++)
			{
				complex <double> comp(0, -2 * MATH_PI*v*y / imgTmp.m_xsize);
				fxv[x*imgTmp.m_xsize + v] += imgTmp.m_data[x*imgTmp.m_xsize + y] * 
					exp(comp);
			}
			fxv[x*imgTmp.m_xsize + v] /= imgTmp.m_xsize;
		} 
	}
	//对矩阵 Fxv 的每一列做 DFT
	for (int v = 0; v < imgTmp.m_xsize; v++)
	{
		for (int u = 0; u < imgTmp.m_ysize; u++)
		{
			fuv[u*imgTmp.m_xsize + v] = complex <double>(0, 0);
			for (int x = 0; x < imgTmp.m_ysize; x++)
			{
				complex <double> comp(0, -2 * MATH_PI*u*x / imgTmp.m_ysize);
				fuv[u*imgTmp.m_xsize + v] += fxv[x*imgTmp.m_xsize + v] * exp(comp);
			}
			fuv[u*imgTmp.m_xsize + v] /= imgTmp.m_ysize;
		} 
	}
	if (FALSE == imgOut1.create(imgIn.m_xsize, imgIn.m_ysize, 3))// DIFT,ERROR
	{
		return FALSE;
	}
	double *data_1 = imgOut1.m_data;
	for(int p=0;p<imgTmp.m_xsize*imgTmp.m_ysize;p++)
	{
		data_1[p] = abs(fuv[p]);//傅立叶谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize] = arg(fuv[p]);//角谱
		data_1[p+imgTmp.m_xsize*imgTmp.m_ysize*2]=norm(fuv[p]);//能量谱
	}
	//对影像上半部分去除条带
	for (int row = 0; row < imgIn.m_ysize / 2; row++)
	{
		for (int col = imgIn.m_xsize / 2; col < imgIn.m_xsize; col++)
		{
			//已知中心辐射线与水平轴的夹角为 78.99，在此设置角度取值在[73.99,83.99]之间
			if (atan2((imgIn.m_ysize / 2.0 - row), (col - imgIn.m_xsize / 2.0)) < (83.99*MATH_PI / 
				180.0) && atan2((imgIn.m_ysize / 2.0 - row), (col - imgIn.m_xsize / 2.0)) > (73.99*MATH_PI / 180.0))
				fuv[row*imgIn.m_ysize + col] = 0;
		} 
	}
	//对影像下半部分去除条带
	for (int row = imgIn.m_ysize / 2; row < imgIn.m_ysize; row++)
	{
		for (int col = 0; col < imgIn.m_xsize / 2; col++)
		{
			//已知中心辐射线与水平轴的夹角为 78.99，在此设置角度取值在[73.99,83.99]之间
			if (atan2((row - imgIn.m_ysize / 2.0), (imgIn.m_xsize / 2.0 - col)) < (83.99*MATH_PI / 
				180.0) && atan2((row - imgIn.m_ysize / 2.0), (imgIn.m_xsize / 2.0 - col)) > (73.99*MATH_PI / 180.0))
				fuv[row*imgIn.m_ysize + col] = 0;
		} 
	}
	//对图像每一行做一维 IDFT
	for (int u = 0; u < imgTmp.m_ysize; u++)
	{
		for (int y = 0; y < imgTmp.m_xsize; y++)
		{
			fuy[u*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int y1 = 0; y1 < imgTmp.m_xsize; y1++)
			{
				complex <double> comp(0, 2 * MATH_PI*y*y1 / imgTmp.m_xsize);
				fuy[u*imgTmp.m_xsize + y] += fuv[u*imgTmp.m_xsize + y1] * exp(comp);
			}
		} 
	}
	//对矩阵 Fvy 的每一列做 IDFT
	for (int y = 0; y < imgTmp.m_xsize; y++)
	{
		for (int x = 0; x < imgTmp.m_ysize; x++)
		{
			fxy[x*imgTmp.m_xsize + y] = complex <double>(0, 0);
			for (int x1 = 0; x1 < imgTmp.m_ysize; x1++)
			{
				complex <double> comp(0, 2 * MATH_PI*x*x1 / imgTmp.m_ysize);
				fxy[x*imgTmp.m_xsize + y] += fuy[x1*imgTmp.m_xsize + y] * exp(comp);
			}
		} 
	}
	if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 2))
	{
		return FALSE;
	}
	double *data = imgOut.m_data;
	for (row = 0; row < imgOut.m_ysize; row++)
	{
		for (col = 0; col < imgOut.m_xsize; col++)
		{
			if ((row + col) % 2 == 1)
			{
				data[row*imgOut.m_xsize + col] = -1 * fxy[row*imgOut.m_xsize + 
					col].real(); 
			}
			else
			{
				data[row*imgOut.m_xsize + col] = fxy[row*imgOut.m_xsize + col].real();
			} 
		} 
	}
	//计算离散傅里叶变换的结果与测试图像之间的误差
	for (int p = 0; p < imgTmp.m_xsize*imgTmp.m_ysize; p++)
	{
		data[p + imgTmp.m_xsize*imgTmp.m_ysize] = abs(imgIn.m_data[p] - data[p]);
	}
	delete fxv;
	delete fuv;
	delete fuy;
	delete fxy;
	return TRUE;
}

BOOL CImageProcessingEx::PCA(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgIn1, CImageDataset &error, int p, int T)
{
	int i,j,k,q,n;
	//定义 CMatrix 矩阵对象
	CMatrix mtx,mtx2,mtx3;
	//定义矩阵大小
	mtx.create(imgIn.m_rastercount,imgIn.m_rastercount);
	mtx2.create(imgIn.m_rastercount,imgIn.m_xsize*imgIn.m_ysize);
	mtx3.create(imgIn.m_rastercount,imgIn.m_xsize*imgIn.m_ysize);
	double **S=new double *[imgIn.m_rastercount];
	for(i=0;i<imgIn.m_rastercount;i++)
	{
		S[i]=new double[imgIn.m_rastercount];
	}
	/* checking image validation */
	if(imgIn.empty())
	{
		return FALSE;
	}
	/* creating the output image */
	if(FALSE == imgIn.duplicate(imgOut))//get new image from input image
	{
		return FALSE;
	}
	if(FALSE == imgIn.duplicate(imgIn1))//get new image from input image
	{
		return FALSE;
	}
	if(FALSE == imgIn.duplicate(error))//get new image from input image
	{
		return FALSE;
	}
	double *data1 = imgIn.m_data;
	double *data2 = imgOut.m_data;
	double *data3 = imgIn1.m_data;
	/*计算协方差矩阵 S*/
	double sum_xi,sum_xj;
	double mean[7];
	n=imgIn.m_ysize*imgIn.m_xsize;
	for(k=0;k<imgIn.m_rastercount;k++)
	{
		mean[k]=0;
		for(int row=0; row<imgIn.m_ysize; row++)
		{
			for(int col=0; col<imgIn.m_xsize; col++)
			{
				mean[k]+=data1[k*n+row*imgIn.m_xsize+col];//
			} 
		}
		mean[k]/=n;
	}
	for(i=0;i<imgIn.m_rastercount;i++)
	{
		for(j=0;j<imgIn.m_rastercount;j++)
		{
			S[i][j]=0;
			for(k=0;k<imgIn.m_xsize*imgIn.m_ysize;k++)
			{
				S[i][j]+=(data1[i*imgIn.m_xsize*imgIn.m_ysize+k]-
					mean[i])*(data1[j*imgIn.m_xsize*imgIn.m_ysize+k]-mean[j]);
			}
			S[i][j]/=(imgIn.m_xsize*imgIn.m_ysize);
		} 
	}
	for(i=0;i<imgIn.m_rastercount;i++)
	{
		for(j=0;j<imgIn.m_rastercount;j++)
		{
			mtx.m_data[i*imgIn.m_rastercount+j]=S[i][j];
		} 
	}
	//特征分解
	float *eigvalues=(float*)malloc(sizeof(float)*mtx.m_rows);
	float *eigvectors=(float*)malloc(sizeof(float)*mtx.m_rows*mtx.m_cols);
	mtx.eig(eigvalues,eigvectors);
	/*从大到小排列特征值并将对应的特征向量也排序*/
	float max,temp;
	int max_i;
	for(i=0;i<imgIn.m_rastercount-1;i++)
	{
		max=eigvalues[i];
		max_i=i;
		for(j=i;j<imgIn.m_rastercount;j++)
		{
			if(eigvalues[j]>max)
			{
				max=eigvalues[j];
				max_i=j;
			} 
		}
		temp=eigvalues[i];
		eigvalues[i]=eigvalues[max_i];
		eigvalues[max_i]=temp;
		for(k=0;k<imgIn.m_rastercount;k++)
		{
			temp=eigvectors[i*imgIn.m_rastercount+k];
			eigvectors[i*imgIn.m_rastercount+k]=eigvectors[max_i*imgIn.m_rastercount+k];
			eigvectors[max_i*imgIn.m_rastercount+k]=temp;
		} 
	}
	double sum;
	/*计算主分量*/
	for(i=0;i<imgIn.m_rastercount;i++)
	{
		for(j=0;j<imgIn.m_xsize*imgIn.m_ysize;j++)
		{
			sum=0;
			for(k=0;k<imgIn.m_rastercount;k++)
			{
				sum+=(eigvectors[i*imgIn.m_rastercount+k]*data1[j+k*imgIn.m_xsize*imgIn.m_ysize]);
			}
			mtx2.m_data[i*imgIn.m_xsize*imgIn.m_ysize+j]=sum;
		} 
	}
	/*提取第 T 个主成分影像*/
	for(i=0;i<imgOut.m_ysize;i++)
	{
		for(j=0;j<imgOut.m_xsize;j++)
		{
			data2[i*imgOut.m_xsize+j]=mtx2.m_data[(T-1)*imgIn.m_xsize*imgIn.m_ysize+i*imgOut.m_xsize+j];
		} 
	}
	/*反变换*/
	for(i=0;i<imgIn.m_rastercount;i++)
	{
		for(j=0;j<imgIn.m_xsize*imgIn.m_ysize;j++)
		{
			for(k=0;k<p;k++)
			{
				mtx3.m_data[i*imgIn.m_xsize*imgIn.m_ysize+j]+=(eigvectors[k*imgIn.m_rastercount+i]*mtx2.m_data[k*imgIn.m_xsize*imgIn.m_ysize+j]);
			} 
		} 
	}
	/*反变换后的矩阵恢复为图像*/
	for(i=0;i<imgOut.m_rastercount;i++)
	{
		for(j=0;j<imgOut.m_xsize*imgOut.m_ysize;j++)
		{
			data3[i*imgOut.m_xsize*imgOut.m_ysize+j]=mtx3.m_data[i*imgOut.m_xsize*imgOut.m_ysize+j];
		} 
	}
	/*计算误差图*/
	double *data_error = error.m_data;
	for(int i=0;i<error.m_rastercount;i++)
	{
		for (int p = 0; p < error.m_xsize*error.m_ysize; p++)
		{
			data_error[p + error.m_xsize*error.m_ysize*i] = abs(imgIn.m_data[p + error.m_xsize*error.m_ysize*i] - imgIn1.m_data [p + error.m_xsize*error.m_ysize*i]);
		}
	}
	//释放 S 数组内存
	for(i=0;i<imgIn.m_rastercount;i++)
	{
		delete[] S[i];
	}
	delete[] S;
	//释放放置特征分解结果的两个内存块
	free(eigvalues);
	free(eigvectors);

	return TRUE; 
}

BOOL CImageProcessingEx::B_I(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgOut1, CImageDataset &error)
{

	int row, col,p,n,min_rgb,s_rgb;
	double fenzi,fenmu,xita,I,H,S;
	/* checking image validation */
	if(imgIn.empty())
	{
		return FALSE;
	}
	/* creating the output image */
	if(FALSE == imgIn.duplicate(imgOut))//get new image from input image
	{
		return FALSE;
	}
	if(FALSE == imgIn.duplicate(imgOut1))//get new image from input image
	{
		return FALSE;
	}
	if(FALSE == imgIn.duplicate(error))//get new image from input image
	{
		return FALSE;
	}
	//RGB到ISH
	double *data = imgIn.m_data;
	n=imgOut.m_ysize*imgOut.m_xsize;
	for(row=0; row<imgOut.m_ysize; row++)
	{
		for(col=0; col<imgOut.m_xsize; col++)
		{
			s_rgb=0;
			for(p=0;p<imgOut.m_rastercount;p++)
			{
				s_rgb=s_rgb+data[p*n+row*imgOut.m_xsize+col];//三通道各像素值累加
			}
			imgOut.m_data[row*imgOut.m_xsize+col]=s_rgb/3;//计算I
			if((data[row*imgOut.m_xsize+col]==data[n+row*imgOut.m_xsize+col])&&(data[n+row*imgOut.m_xsize+col]==data[2*n+row*imgOut.m_xsize+col]))
				imgOut.m_data[n+row*imgOut.m_xsize+col]=0;//特殊情况
			else
			{
				fenzi=(data[row*imgOut.m_xsize+col]-data[n+row*imgOut.m_xsize+col]+data[row*imgOut.m_xsize+col]-data[2*n+row*imgOut.m_xsize+col])/2;
				fenmu=pow((pow(data[row*imgOut.m_xsize+col]-data[n+row*imgOut.m_xsize+col],2)+(data[row*imgOut.m_xsize+col]-data[2*n+row*imgOut.m_xsize+col])*(data[n+row*imgOut.m_xsize+col]-data[2*n+row*imgOut.m_xsize+col])),0.5);
				xita=acos(fenzi/fenmu);
				if(data[2*n+row*imgOut.m_xsize+col]>data[n+row*imgOut.m_xsize+col])//计算 H
					imgOut.m_data[n+row*imgOut.m_xsize+col]=360-xita*180/pi;
				else
					imgOut.m_data[n+row*imgOut.m_xsize+col]=xita*180/pi;
			}
			if(s_rgb==0)
				imgOut.m_data[2*n+row*imgOut.m_xsize+col]=0;
			else
			{ /*三者比较找到最小值*/
				min_rgb=data[row*imgOut.m_xsize+col];
				if (min_rgb>data[n+row*imgOut.m_xsize+col])
					min_rgb=data[n+row*imgOut.m_xsize+col];
				if (min_rgb>data[2*n+row*imgOut.m_xsize+col])
					min_rgb=data[2*n+row*imgOut.m_xsize+col];
				imgOut.m_data[2*n+row*imgOut.m_xsize+col]=1-3*min_rgb/s_rgb;//计算 S 
			} 
		} 
	}

	//IHS变换到RGB
	double *data_1 = imgOut.m_data;
	for(row=0; row<imgOut1.m_ysize; row++)
	{
		for(col=0; col<imgOut1.m_xsize; col++)
		{
			I=data_1[row*imgOut1.m_xsize+col];//幅值给 IHS 简化变量
			H=data_1[n+row*imgOut1.m_xsize+col]*pi/180;
			S=data_1[2*n+row*imgOut1.m_xsize+col];
			if((H>=0)&&(H<(120*pi/180)))//if 条件语句，控制代入公式
			{
				imgOut1.m_data[row*imgOut1.m_xsize+col]=I*(1+(S*cos(H))/(cos(60*pi/180-H)));//计算 R
				imgOut1.m_data[2*n+row*imgOut1.m_xsize+col]=I*(1-S);//计 算 B
				imgOut1.m_data[n+row*imgOut1.m_xsize+col]=3*I-(imgOut1.m_data[row*imgOut1.m_xsize+col]+imgOut1.m_data[2*n+row*imgOut1.m_xsize+col]);//计算 G 
			}
			else if(H<(240*pi/180))
			{
				imgOut1.m_data[row*imgOut1.m_xsize+col]=I*(1-S);//计算 R
				imgOut1.m_data[n+row*imgOut1.m_xsize+col]=I*(1+((S*cos(H-(pi*120/180)))/(cos(pi-H))));//计算 G
				imgOut1.m_data[2*n+row*imgOut1.m_xsize+col]=3*I-(imgOut1.m_data[row*imgOut1.m_xsize+col]+imgOut1.m_data[n+row*imgOut1.m_xsize+col]);//计算 B 
			}
			else if(H<2*pi)
			{
				imgOut1.m_data[n+row*imgOut1.m_xsize+col]=I*(1-S);//计算G
				imgOut1.m_data[2*n+row*imgOut1.m_xsize+col]=I*(1+((S*cos(H-(pi*240/180)))/(cos((300*pi/180)-H))));//计算 B
				imgOut1.m_data[row*imgOut1.m_xsize+col]=3*I-(imgOut1.m_data[n+row*imgOut1.m_xsize+col]+imgOut1.m_data[2*n+row*imgOut1.m_xsize+col]);//计算 R 
			} 
		}
	}

	//计算误差图
	double *data_error = error.m_data;
	for(int i=0;i<error.m_rastercount;i++)
	{
		for (int p = 0; p < error.m_xsize*error.m_ysize; p++)
		{
			data_error[p + error.m_xsize*error.m_ysize*i] = abs(imgIn.m_data[p + error.m_xsize*error.m_ysize*i] - imgOut1.m_data [p + error.m_xsize*error.m_ysize*i]);
		}
	}
	return TRUE;
}

//BOOL CImageProcessingEx::I_B(CImageDataset &imgIn, CImageDataset &imgOut)
//{
//	int row, col,p,n;
//	double I,H,S;
//	/* checking image validation */
//	if(imgIn.empty())
//	{
//		return FALSE;
//	}/* creating the output image */
//	if(FALSE == imgIn.duplicate(imgOut))//get new image from input image
//	{
//		return FALSE;
//	}
//	double *data = imgIn.m_data;
//	n=imgOut.m_ysize*imgOut.m_xsize;
//	for(row=0; row<imgOut.m_ysize; row++)
//	{
//		for(col=0; col<imgOut.m_xsize; col++)
//		{
//			I=data[row*imgOut.m_xsize+col];//幅值给 IHS 简化变量
//			H=data[n+row*imgOut.m_xsize+col]*pi/180;
//			S=data[2*n+row*imgOut.m_xsize+col];
//			if((H>=0)&&(H<(120*pi/180)))//if 条件语句，控制代入公式
//			{
//				imgOut.m_data[row*imgOut.m_xsize+col]=I*(1+(S*cos(H))/(cos(60*pi/180-H)));//计算 R
//				imgOut.m_data[2*n+row*imgOut.m_xsize+col]=I*(1-S);//计 算 B
//				imgOut.m_data[n+row*imgOut.m_xsize+col]=3*I-(imgOut.m_data[row*imgOut.m_xsize+col]+imgOut.m_data[2*n+row*imgOut.m_xsize+col]);//计算 G 
//			}
//			else if(H<(240*pi/180))
//			{
//				imgOut.m_data[row*imgOut.m_xsize+col]=I*(1-S);//计算 R
//				imgOut.m_data[n+row*imgOut.m_xsize+col]=I*(1+((S*cos(H-(pi*120/180)))/(cos(pi-H))));//计算 G
//				imgOut.m_data[2*n+row*imgOut.m_xsize+col]=3*I-(imgOut.m_data[row*imgOut.m_xsize+col]+imgOut.m_data[n+row*imgOut.m_xsize+col]);//计算 B 
//			}
//			else if(H<=2*pi)
//			{
//				imgOut.m_data[n+row*imgOut.m_xsize+col]=I*(1-S);//计算G
//				imgOut.m_data[2*n+row*imgOut.m_xsize+col]=I*(1+((S*cos(H-(pi*240/180)))/(cos((300*pi/180)-H))));//计算 B
//				imgOut.m_data[row*imgOut.m_xsize+col]=3*I-(imgOut.m_data[n+row*imgOut.m_xsize+col]+imgOut.m_data[2*n+row*imgOut.m_xsize+col]);//计算 R 
//			} 
//		}
//	}
//}