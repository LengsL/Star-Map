
// starmap2Dlg.h : header file
//

#pragma once


// Cstarmap2Dlg dialog
class Cstarmap2Dlg : public CDialogEx
{
// Construction
public:
	Cstarmap2Dlg(CWnd* pParent = nullptr);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_STARMAP2_DIALOG };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	virtual BOOL OnInitDialog();

	DECLARE_MESSAGE_MAP()
private:
	HICON m_hIcon;
    CString m_gotoResultPath;
    ULONGLONG m_lastGotoSequence = 0;

    bool ReadCoordinate(
        int controlId,
        LPCTSTR label,
        double minimum,
        double maximum,
        bool includeMaximum,
        double& value) const;

    bool LaunchPythonStarMap(
        double altitude,
        double azimuth,
        double longitude,
        double latitude,
        double height);

    void ReadGotoResult();

public:
    afx_msg void OnBnClickedButtonDisplay();
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnDestroy();
	
};
