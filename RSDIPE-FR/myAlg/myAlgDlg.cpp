
// myAlgDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "myAlg.h"
#include "myAlgDlg.h"
#include "afxdialogex.h"
#include "math.h"

#include "ImageProcessingEx.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CmyAlgDlg 对话框




CmyAlgDlg::CmyAlgDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CmyAlgDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
	m_strimginput = _T("");
	m_strimgoutput = _T("");
	//  m_size = 0;
	m_size = 0;
	m_N = 0;
	m_R = 0.0f;
	m_S = 0.0f;
	m_P = 0;
	m_T = 0.0f;
	m_BD0 = 0.0f;
	m_BN = 0;
	m_GD0 = 0.0f;
	m_LD0 = 0.0f;
}

void CmyAlgDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, m_strimginput);
	DDX_Text(pDX, IDC_EDIT2, m_strimgoutput);
	//  DDX_Text(pDX, IDC_EDIT_MID, m_size);
	DDX_Text(pDX, IDC_EDIT_MID, m_size);
	DDX_Text(pDX, IDC_EDIT_N, m_N);
	DDX_Text(pDX, IDC_EDIT_R, m_R);
	DDX_Text(pDX, IDC_EDIT_S, m_S);
	DDX_Text(pDX, IDC_EDIT_P, m_P);
	DDX_Text(pDX, IDC_EDIT_T, m_T);
	DDX_Text(pDX, IDC_EDIT_BD0, m_BD0);
	DDX_Text(pDX, IDC_EDIT_BN, m_BN);
	DDX_Text(pDX, IDC_EDIT_GD0, m_GD0);
	DDX_Text(pDX, IDC_EDIT_LD0, m_LD0);
}

BEGIN_MESSAGE_MAP(CmyAlgDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_READ, &CmyAlgDlg::OnClickedButtonRead)
	ON_BN_CLICKED(IDC_BUTTON_SAVE, &CmyAlgDlg::OnClickedButtonSave)
	ON_BN_CLICKED(IDC_BUT_HISTEQ, &CmyAlgDlg::OnClickedButHisteq)
	ON_BN_CLICKED(IDC_BUT_MATCH, &CmyAlgDlg::OnClickedButMatch)
	ON_BN_CLICKED(IDC_BUTTON_MID, &CmyAlgDlg::OnClickedButtonMid)
	ON_EN_CHANGE(IDC_EDIT_MID, &CmyAlgDlg::OnEnChangeEditMid)
	ON_BN_CLICKED(IDC_BUTTON_BIL, &CmyAlgDlg::OnClickedButtonBil)
	//ON_EN_CHANGE(IDC_EDIT_R, &CmyAlgDlg::OnEnChangeEd/*itR)
	//ON_EN_CHANGE(IDC_EDIT_S, &CmyAlgDlg::OnEnChangeEditS)*/
	ON_BN_CLICKED(IDC_BUTTON_LAP, &CmyAlgDlg::OnClickedButtonLap)
//	ON_BN_CLICKED(IDC_BUTTON_DFT, &CmyAlgDlg::OnClickedButtonDft)

	ON_BN_CLICKED(IDC_BUTTON_Creat, &CmyAlgDlg::OnClickedButtonCreat)

	//ON_BN_CLICKED(IDC_BUTTON_IDFT, &CmyAlgDlg::OnBnClickedButtonIdft)
	ON_BN_CLICKED(IDC_BUTTON_LG, &CmyAlgDlg::OnBnClickedButtonLg)
	ON_BN_CLICKED(IDC_BUTTON_BTWG, &CmyAlgDlg::OnBnClickedButtonBtwg)

	ON_BN_CLICKED(IDC_BUTTON_GG, &CmyAlgDlg::OnBnClickedButtonGg)
	ON_BN_CLICKED(IDC_BUTTON_QTD, &CmyAlgDlg::OnBnClickedButtonQtd)
	ON_BN_CLICKED(IDC_BUTTON_PCA, &CmyAlgDlg::OnBnClickedButtonPca)
	ON_BN_CLICKED(IDC_BUTTON_R_I, &CmyAlgDlg::OnBnClickedButtonRI)
	//ON_BN_CLICKED(IDC_BUTTON_I_R, &CmyAlgDlg::OnBnClickedButtonIR)

	ON_BN_CLICKED(IDC_BUTTON_DFT, &CmyAlgDlg::OnClickedButtonDft)
END_MESSAGE_MAP()


// CmyAlgDlg 消息处理程序

BOOL CmyAlgDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码

	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CmyAlgDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CmyAlgDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CmyAlgDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CmyAlgDlg::OnClickedButtonRead()
{
	CFileDialog dlg(TRUE);///TRUE 为 OPEN 对话框，FALSE 为 SAVE AS 对话框 
  if(dlg.DoModal()==IDOK)  
  {  
    if(FALSE == CImageIO::read(imgIn, dlg.GetPathName()))  
    {  
      AfxMessageBox("读取图像失败！");  
    }  
    else  
    {  
      m_strimginput = dlg.GetPathName();  
      UpdateData(FALSE);  
      if(imgIn.m_rastercount>=3)  
        CImageDisplay::show(imgIn, this, dlg.GetFileName(), 1, 2, 3, 0);  
      else  
        CImageDisplay::show(imgIn, this, dlg.GetFileName(), 1, 1, 1, 0);  
	}

  }
}


void CmyAlgDlg::OnClickedButtonSave()
{
	CFileDialog dlg(FALSE);///TRUE 为 OPEN 对话框，FALSE 为 SAVE AS 对话框 
	if(dlg.DoModal()==IDOK) 
	{ 
		if(FALSE == CImageIO::write(imgOut, dlg.GetPathName())) 
		{ 
			AfxMessageBox("保存图像失败！"); 
		} 
		else 
		{ 
			m_strimgoutput = dlg.GetPathName(); 
			UpdateData(FALSE); 
			if(imgOut.m_rastercount>=3) 
				CImageDisplay::show(imgOut, this, dlg.GetFileName(), 1, 2, 3, 0); 
			else 
				CImageDisplay::show(imgOut, this, dlg.GetFileName(), 1, 1, 1, 0); 
		}
	}
}



void CmyAlgDlg::OnClickedButHisteq()
{
  int bSuccess = CImageProcessingEx::histeq(imgIn, imgOut);  
  if(!bSuccess)  
  {  
    AfxMessageBox("处理图像失败！");  
    return;  
  }  
  if(imgOut.m_rastercount>=3)  
    CImageDisplay::show(imgOut, this, "直方图均衡化", 1, 2, 3, 0);  
  else  
    CImageDisplay::show(imgOut, this, "直方图均衡化", 1, 1, 1, 0);  
} 



void CmyAlgDlg::OnClickedButMatch()
{
	CImageDataset imgRefer;
	CFileDialog dlg(TRUE);///TRUE 为 OPEN 对话框，FALSE 为 SAVE AS 对话框 
	if(dlg.DoModal()==IDOK) 
	{ 
		if(FALSE == CImageIO::read(imgRefer, dlg.GetPathName())) 
		{ 
			AfxMessageBox("读取图像失败！"); 
		} 
		else 
		{ 
			//m_strimginput = dlg.GetPathName();
			UpdateData(FALSE); 
			if(imgIn.m_rastercount>=3) 
				CImageDisplay::show(imgRefer, this, "参考图像", 1, 2, 3, 0); 
			else 
				CImageDisplay::show(imgRefer, this, "参考图像", 1, 1, 1, 0); 
			int bSuccess = CImageProcessingEx::match(imgIn,imgRefer, imgOut); 
			if(!bSuccess) 
			{ 
				AfxMessageBox("处理图像失败！"); 
				return; 
			} 
			if(imgOut.m_rastercount>=3) 
				CImageDisplay::show(imgOut, this, "直方图匹配", 1, 2, 3, 0); 
			else 
				CImageDisplay::show(imgOut, this, "直方图匹配", 1, 1, 1, 0); 
		} 
	} 

}


void CmyAlgDlg::OnClickedButtonMid()
{
	UpdateData(TRUE);
	int m=m_size;
	if (m<3 || (m%2)==0)
	{
		AfxMessageBox("错误！请重新输入滤波器窗口大小！");
		return;
	}

	int bSuccess = CImageProcessingEx::median(m,imgIn, imgOut);

	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	if(imgOut.m_rastercount>=3)
		CImageDisplay::show(imgOut, this, "中值滤波", 1, 2, 3, 0);
	else
		CImageDisplay::show(imgOut, this, "中值滤波", 1, 1, 1, 0);
}


void CmyAlgDlg::OnEnChangeEditMid()
{
	// TODO:  如果该控件是 RICHEDIT 控件，它将不
	// 发送此通知，除非重写 CDialogEx::OnInitDialog()
	// 函数并调用 CRichEditCtrl().SetEventMask()，
	// 同时将 ENM_CHANGE 标志“或”运算到掩码中。

	// TODO:  在此添加控件通知处理程序代码
}


void CmyAlgDlg::OnClickedButtonBil()
{
	UpdateData(TRUE);//更新变量
	int N1=m_N;
	if (m_N<3 || (m_N%2)==0)
	{
		AfxMessageBox("参数错误！请重新输入滤波器窗口大小！");
		return;
	}
	double s1=m_S;
	double r1=m_R;

	int bSuccess = CImageProcessingEx::bilateral(N1,s1,r1,imgIn, imgOut);//将在页面上输入的参数作为函数参数
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	if(imgOut.m_rastercount>=3)
		CImageDisplay::show(imgOut, this, "双边滤波", 1, 2, 3, 0);
	else
		CImageDisplay::show(imgOut, this, "双边滤波", 1, 1, 1, 0);
}


void CmyAlgDlg::OnClickedButtonLap()
{
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::laplacian(imgIn, imgOut);
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	if(imgOut.m_rastercount>=3)
		CImageDisplay::show(imgOut, this, "锐化处理", 1, 2, 3, 0);
	else
		CImageDisplay::show(imgOut, this, "锐化处理", 1, 1, 1, 0);
}


void CmyAlgDlg::OnClickedButtonCreat()
{
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::Creat(imgOut);
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
	}
	CImageDisplay::show(imgOut, this, "测试图像", 1, 1, 1, 0);
}


void CmyAlgDlg::OnBnClickedButtonLg()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);//更新变量
	int bSuccess = CImageProcessingEx::LG(imgIn, imgOut, imgOut1, m_LD0);//使用页面上输入的变量参数
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT结果", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "对数拉伸", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "角谱", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "能量谱", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "理想高通滤波器", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonBtwg()
{
	UpdateData(TRUE);//更新变量
	int bSuccess = CImageProcessingEx::BTWG(imgIn, imgOut, imgOut1, m_BD0, m_BN);//使用页面上输入的变量参数
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT结果", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "对数拉伸", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "角谱", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "能量谱", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "巴特沃斯高通滤波器", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonGg()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);//更新变量
	int bSuccess = CImageProcessingEx::GG(imgIn, imgOut, imgOut1, m_GD0);//使用页面上输入的变量参数
	if(!bSuccess) {
		AfxMessageBox("处理图像失败！");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT结果", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "对数拉伸", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "角谱", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "能量谱", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "高斯高通滤波器", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonQtd()
{
	// TODO: 在此添加控件通知处理程序代码
	int bSuccess = CImageProcessingEx::QTD(imgIn, imgOut, imgOut1);
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT结果", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "对数拉伸", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "角谱", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "能量谱", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "去除条带", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonPca()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::PCA(imgIn, imgOut,imgIn1,error,m_P,m_T);
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	CImageDisplay::show(imgOut, this, "第T个主分量图", 1, 1, 1, 0);
	CImageDisplay::show(imgIn1, this, "pca 反变换", 1, 2, 3, 0);
	/*计算重建后的遥感影像与原始遥感影像之间的误差*/
	CImageDisplay::show(error, this, "pca 误差图", 1, 2, 3, 1);

	double errors=0;
	double *data=error.m_data;
	//double *dataOut=imgOut.m_data;
	for(int row=0; row<imgIn.m_ysize; row++)
	{
		for(int col=0; col<imgIn.m_xsize; col++)
		{
			errors+=(data[row*imgIn.m_xsize+col])*(data[row*imgIn.m_xsize+col]);
		}
	}
	errors=sqrt(errors);

	CString str;
	str.Format("均方根误差为：%f", errors);
	AfxMessageBox(str);
}


void CmyAlgDlg::OnBnClickedButtonRI()
{
	// TODO: 在此添加控件通知处理程序代码
	int bSuccess = CImageProcessingEx::B_I(imgIn, imgOut, imgOut1, error);
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	if(imgOut.m_rastercount>=3)
	{
		CImageDisplay::show(imgOut, this, "RGB-IHS 未拉伸", 1, 2, 3, 0);
		CImageDisplay::show(imgOut, this, "RGB-IHS 拉伸", 1, 2, 3, 1);
	}
	else
		CImageDisplay::show(imgOut, this, "RGB-IHS", 1, 1, 1, 0);

	if(imgOut1.m_rastercount>=3)
	{
		CImageDisplay::show(imgOut1, this, "IHS-RGB", 1, 2, 3, 0);
	}
	else
		CImageDisplay::show(imgOut1, this, "IHS-RGB", 1, 1, 1, 0);

	/*计算重建后的遥感影像与原始遥感影像之间的误差*/
	CImageDisplay::show(error, this, "误差图", 1, 2, 3, 1);

	double errors=0;
	double *data=error.m_data;
	//double *dataOut=imgOut.m_data;
	for(int row=0; row<imgIn.m_ysize; row++)
	{
		for(int col=0; col<imgIn.m_xsize; col++)
		{
			errors+=(data[row*imgIn.m_xsize+col])*(data[row*imgIn.m_xsize+col]);
		}
	}
	errors=sqrt(errors);

	CString str;
	str.Format("均方根误差为：%f", errors);
	AfxMessageBox(str);

}


void CmyAlgDlg::OnClickedButtonDft()
{
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::DFT(imgIn, imgOut, imgOut1, error);
	if(!bSuccess)
	{
		AfxMessageBox("处理图像失败！");
		return;
	}
	CImageDisplay::show(imgOut, this, "二维傅里叶反变换", 1, 1, 1, 0);
	for(int i=0;i<imgOut.m_ysize;i++)
	{
		for(int j=0;j<imgOut.m_xsize;j++)
		{
			imgOut.m_data[i*imgOut.m_xsize+j]=fabs(imgOut.m_data[i*imgOut.m_xsize+j]-imgIn.m_data[i*imgOut.m_xsize+j]);
		} 
	}
	
	CImageDisplay::show(error, this, "误差", 1, 1, 1, 1);

	double errors=0;
	double *data=error.m_data;
	//double *dataOut=imgOut.m_data;
	for(int row=0; row<imgIn.m_ysize; row++)
	{
		for(int col=0; col<imgIn.m_xsize; col++)
		{
			errors+=(data[row*imgIn.m_xsize+col])*(data[row*imgIn.m_xsize+col]);
		}
	}
	errors=sqrt(errors);

	CString str;
	str.Format("均方根误差为：%f", errors);
	AfxMessageBox(str);

	CImageDisplay::show(imgOut1, this, "DFT结果", 1, 1, 1, 0);
	CImageDisplay::show(imgOut1, this, "线性变换（0-255）", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "对数拉伸", 1, 1, 1, 3);
	//CImageDisplay::show(imgOut1, this, "能量谱", 3, 3, 3, 3);
}
