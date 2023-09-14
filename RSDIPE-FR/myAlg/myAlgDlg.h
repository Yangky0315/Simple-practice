
// myAlgDlg.h : 头文件
//

#pragma once
#include "RSDIPLib.h"


// CmyAlgDlg 对话框
class CmyAlgDlg : public CDialogEx
{
// 构造
public:
	CmyAlgDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_MYALG_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
private:
	afx_msg void OnClickedButtonRead();
	CString m_strimginput;

	CImageDataset imgIn;
	CImageDataset imgRefer;
	CImageDataset imgIn1;
	CImageDataset imgOut;
	CImageDataset imgOut1;
	CImageDataset error;
	CImageDataset imgNew;
	afx_msg void OnClickedButtonSave();
	CString m_strimgoutput;
	afx_msg void OnClickedButHisteq();
	afx_msg void OnClickedButMatch();
	afx_msg void OnBnClickedButton1();
	afx_msg void OnClickedButtonMid();

	int m_size;
	afx_msg void OnEnChangeEditMid();
	afx_msg void OnClickedButtonBil();
	int m_N;
	float m_R;
	float m_S;
	//afx_msg void OnEnChangeEditR();/*
	//afx_msg void OnEnChangeEditS();*/
	afx_msg void OnClickedButtonLap();
//	afx_msg void OnClickedButtonDft();
	afx_msg void OnClickedButtonNew();
	afx_msg void OnClickedButtonCreat();
	afx_msg void OnBnClickedButtonG();
	//afx_msg void OnBnClickedButtonIdft();
	afx_msg void OnBnClickedButtonLg();
	afx_msg void OnBnClickedButtonBtwg();
	afx_msg void OnBnClickedButton3();
	afx_msg void OnBnClickedButtonGg();
	afx_msg void OnBnClickedButtonQtd();
	int m_P;
	afx_msg void OnBnClickedButtonPca();
	float m_T;
	afx_msg void OnBnClickedButtonRI();
	//afx_msg void OnBnClickedButtonIR();
	afx_msg void OnEnChangeEditD0();
	float m_BD0;
	int m_BN;
	float m_GD0;
	float m_LD0;
public:
	afx_msg void OnClickedButtonDft();
};
