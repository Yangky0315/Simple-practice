
// myAlgDlg.cpp : ʵ���ļ�
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


// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// �Ի�������
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

// ʵ��
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


// CmyAlgDlg �Ի���




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


// CmyAlgDlg ��Ϣ�������

BOOL CmyAlgDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// ��������...���˵�����ӵ�ϵͳ�˵��С�

	// IDM_ABOUTBOX ������ϵͳ���Χ�ڡ�
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

	// ���ô˶Ի����ͼ�ꡣ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, FALSE);		// ����Сͼ��

	// TODO: �ڴ���Ӷ���ĳ�ʼ������

	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
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

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void CmyAlgDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR CmyAlgDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CmyAlgDlg::OnClickedButtonRead()
{
	CFileDialog dlg(TRUE);///TRUE Ϊ OPEN �Ի���FALSE Ϊ SAVE AS �Ի��� 
  if(dlg.DoModal()==IDOK)  
  {  
    if(FALSE == CImageIO::read(imgIn, dlg.GetPathName()))  
    {  
      AfxMessageBox("��ȡͼ��ʧ�ܣ�");  
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
	CFileDialog dlg(FALSE);///TRUE Ϊ OPEN �Ի���FALSE Ϊ SAVE AS �Ի��� 
	if(dlg.DoModal()==IDOK) 
	{ 
		if(FALSE == CImageIO::write(imgOut, dlg.GetPathName())) 
		{ 
			AfxMessageBox("����ͼ��ʧ�ܣ�"); 
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
    AfxMessageBox("����ͼ��ʧ�ܣ�");  
    return;  
  }  
  if(imgOut.m_rastercount>=3)  
    CImageDisplay::show(imgOut, this, "ֱ��ͼ���⻯", 1, 2, 3, 0);  
  else  
    CImageDisplay::show(imgOut, this, "ֱ��ͼ���⻯", 1, 1, 1, 0);  
} 



void CmyAlgDlg::OnClickedButMatch()
{
	CImageDataset imgRefer;
	CFileDialog dlg(TRUE);///TRUE Ϊ OPEN �Ի���FALSE Ϊ SAVE AS �Ի��� 
	if(dlg.DoModal()==IDOK) 
	{ 
		if(FALSE == CImageIO::read(imgRefer, dlg.GetPathName())) 
		{ 
			AfxMessageBox("��ȡͼ��ʧ�ܣ�"); 
		} 
		else 
		{ 
			//m_strimginput = dlg.GetPathName();
			UpdateData(FALSE); 
			if(imgIn.m_rastercount>=3) 
				CImageDisplay::show(imgRefer, this, "�ο�ͼ��", 1, 2, 3, 0); 
			else 
				CImageDisplay::show(imgRefer, this, "�ο�ͼ��", 1, 1, 1, 0); 
			int bSuccess = CImageProcessingEx::match(imgIn,imgRefer, imgOut); 
			if(!bSuccess) 
			{ 
				AfxMessageBox("����ͼ��ʧ�ܣ�"); 
				return; 
			} 
			if(imgOut.m_rastercount>=3) 
				CImageDisplay::show(imgOut, this, "ֱ��ͼƥ��", 1, 2, 3, 0); 
			else 
				CImageDisplay::show(imgOut, this, "ֱ��ͼƥ��", 1, 1, 1, 0); 
		} 
	} 

}


void CmyAlgDlg::OnClickedButtonMid()
{
	UpdateData(TRUE);
	int m=m_size;
	if (m<3 || (m%2)==0)
	{
		AfxMessageBox("���������������˲������ڴ�С��");
		return;
	}

	int bSuccess = CImageProcessingEx::median(m,imgIn, imgOut);

	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	if(imgOut.m_rastercount>=3)
		CImageDisplay::show(imgOut, this, "��ֵ�˲�", 1, 2, 3, 0);
	else
		CImageDisplay::show(imgOut, this, "��ֵ�˲�", 1, 1, 1, 0);
}


void CmyAlgDlg::OnEnChangeEditMid()
{
	// TODO:  ����ÿؼ��� RICHEDIT �ؼ���������
	// ���ʹ�֪ͨ��������д CDialogEx::OnInitDialog()
	// ���������� CRichEditCtrl().SetEventMask()��
	// ͬʱ�� ENM_CHANGE ��־�������㵽�����С�

	// TODO:  �ڴ���ӿؼ�֪ͨ����������
}


void CmyAlgDlg::OnClickedButtonBil()
{
	UpdateData(TRUE);//���±���
	int N1=m_N;
	if (m_N<3 || (m_N%2)==0)
	{
		AfxMessageBox("�������������������˲������ڴ�С��");
		return;
	}
	double s1=m_S;
	double r1=m_R;

	int bSuccess = CImageProcessingEx::bilateral(N1,s1,r1,imgIn, imgOut);//����ҳ��������Ĳ�����Ϊ��������
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	if(imgOut.m_rastercount>=3)
		CImageDisplay::show(imgOut, this, "˫���˲�", 1, 2, 3, 0);
	else
		CImageDisplay::show(imgOut, this, "˫���˲�", 1, 1, 1, 0);
}


void CmyAlgDlg::OnClickedButtonLap()
{
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::laplacian(imgIn, imgOut);
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	if(imgOut.m_rastercount>=3)
		CImageDisplay::show(imgOut, this, "�񻯴���", 1, 2, 3, 0);
	else
		CImageDisplay::show(imgOut, this, "�񻯴���", 1, 1, 1, 0);
}


void CmyAlgDlg::OnClickedButtonCreat()
{
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::Creat(imgOut);
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
	}
	CImageDisplay::show(imgOut, this, "����ͼ��", 1, 1, 1, 0);
}


void CmyAlgDlg::OnBnClickedButtonLg()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(TRUE);//���±���
	int bSuccess = CImageProcessingEx::LG(imgIn, imgOut, imgOut1, m_LD0);//ʹ��ҳ��������ı�������
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT���", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "��������", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "����", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "������", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "�����ͨ�˲���", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonBtwg()
{
	UpdateData(TRUE);//���±���
	int bSuccess = CImageProcessingEx::BTWG(imgIn, imgOut, imgOut1, m_BD0, m_BN);//ʹ��ҳ��������ı�������
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT���", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "��������", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "����", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "������", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "������˹��ͨ�˲���", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonGg()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(TRUE);//���±���
	int bSuccess = CImageProcessingEx::GG(imgIn, imgOut, imgOut1, m_GD0);//ʹ��ҳ��������ı�������
	if(!bSuccess) {
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT���", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "��������", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "����", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "������", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "��˹��ͨ�˲���", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonQtd()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	int bSuccess = CImageProcessingEx::QTD(imgIn, imgOut, imgOut1);
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	CImageDisplay::show(imgOut1, this, "DFT���", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "��������", 1, 1, 1, 3);
	CImageDisplay::show(imgOut1, this, "����", 2, 2, 2, 1);
	CImageDisplay::show(imgOut1, this, "������", 3, 3, 3, 3);

	CImageDisplay::show(imgOut, this, "ȥ������", 1, 1, 1, 1);
}


void CmyAlgDlg::OnBnClickedButtonPca()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::PCA(imgIn, imgOut,imgIn1,error,m_P,m_T);
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	CImageDisplay::show(imgOut, this, "��T��������ͼ", 1, 1, 1, 0);
	CImageDisplay::show(imgIn1, this, "pca ���任", 1, 2, 3, 0);
	/*�����ؽ����ң��Ӱ����ԭʼң��Ӱ��֮������*/
	CImageDisplay::show(error, this, "pca ���ͼ", 1, 2, 3, 1);

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
	str.Format("���������Ϊ��%f", errors);
	AfxMessageBox(str);
}


void CmyAlgDlg::OnBnClickedButtonRI()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	int bSuccess = CImageProcessingEx::B_I(imgIn, imgOut, imgOut1, error);
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	if(imgOut.m_rastercount>=3)
	{
		CImageDisplay::show(imgOut, this, "RGB-IHS δ����", 1, 2, 3, 0);
		CImageDisplay::show(imgOut, this, "RGB-IHS ����", 1, 2, 3, 1);
	}
	else
		CImageDisplay::show(imgOut, this, "RGB-IHS", 1, 1, 1, 0);

	if(imgOut1.m_rastercount>=3)
	{
		CImageDisplay::show(imgOut1, this, "IHS-RGB", 1, 2, 3, 0);
	}
	else
		CImageDisplay::show(imgOut1, this, "IHS-RGB", 1, 1, 1, 0);

	/*�����ؽ����ң��Ӱ����ԭʼң��Ӱ��֮������*/
	CImageDisplay::show(error, this, "���ͼ", 1, 2, 3, 1);

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
	str.Format("���������Ϊ��%f", errors);
	AfxMessageBox(str);

}


void CmyAlgDlg::OnClickedButtonDft()
{
	UpdateData(TRUE);
	int bSuccess = CImageProcessingEx::DFT(imgIn, imgOut, imgOut1, error);
	if(!bSuccess)
	{
		AfxMessageBox("����ͼ��ʧ�ܣ�");
		return;
	}
	CImageDisplay::show(imgOut, this, "��ά����Ҷ���任", 1, 1, 1, 0);
	for(int i=0;i<imgOut.m_ysize;i++)
	{
		for(int j=0;j<imgOut.m_xsize;j++)
		{
			imgOut.m_data[i*imgOut.m_xsize+j]=fabs(imgOut.m_data[i*imgOut.m_xsize+j]-imgIn.m_data[i*imgOut.m_xsize+j]);
		} 
	}
	
	CImageDisplay::show(error, this, "���", 1, 1, 1, 1);

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
	str.Format("���������Ϊ��%f", errors);
	AfxMessageBox(str);

	CImageDisplay::show(imgOut1, this, "DFT���", 1, 1, 1, 0);
	CImageDisplay::show(imgOut1, this, "���Ա任��0-255��", 1, 1, 1, 1);
	CImageDisplay::show(imgOut1, this, "��������", 1, 1, 1, 3);
	//CImageDisplay::show(imgOut1, this, "������", 3, 3, 3, 3);
}
