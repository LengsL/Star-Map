// starmap2Dlg.cpp : implementation file
//

#include "pch.h"
#include "framework.h"
#include "starmap2.h"
#include "starmap2Dlg.h"
#include "afxdialogex.h"

#include <atlconv.h>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <Shellapi.h>
#include <tchar.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

namespace
{
    CString GetAppRoot()
    {
        TCHAR buffer[MAX_PATH] = { 0 };
        GetModuleFileName(nullptr, buffer, MAX_PATH);

        CString path(buffer);
        int lastBackslashIndex = path.ReverseFind(_T('\\'));

        if (lastBackslashIndex != -1)
            path = path.Left(lastBackslashIndex);

        return path;
    }
    CString GetStarMapRoot()
    {
        // Support both MFC\x64\Debug and MFC\Debug output directories.
        const CString candidates[] = {
            GetAppRoot() + _T("\\..\\..\\..\\python"),
            GetAppRoot() + _T("\\..\\..\\python"),
            GetAppRoot() + _T("\\..\\python")
        };

        for (const CString& candidate : candidates)
        {
            TCHAR buffer[MAX_PATH] = { 0 };
            GetFullPathName(candidate, MAX_PATH, buffer, nullptr);
            CString root(buffer);
            if (GetFileAttributes(root + _T("\\main.py")) != INVALID_FILE_ATTRIBUTES)
                return root;
        }

        // Return the standard x64 Debug layout for a useful error message.
        TCHAR buffer[MAX_PATH] = { 0 };
        GetFullPathName(candidates[0], MAX_PATH, buffer, nullptr);
        return CString(buffer);
    }

    CString GetMainPyPath()
    {
        return GetStarMapRoot() + _T("\\main.py");
    }

    CString GetPythonExePath()
    {
        return GetStarMapRoot() + _T("\\.venv\\Scripts\\python.exe");
    }

    bool ParseDoubleStrict(
        const CString& source,
        double& value)
    {
        CString text(source);
        text.Trim();

        if (text.IsEmpty())
            return false;

        errno = 0;

        TCHAR* endPtr = nullptr;
        const TCHAR* beginPtr = text.GetString();

        value = _tcstod(beginPtr, &endPtr);

        while (endPtr != nullptr && _istspace(*endPtr))
            ++endPtr;

        return errno != ERANGE &&
               endPtr != beginPtr &&
               endPtr != nullptr &&
               *endPtr == _T('\0');
    }
}


// Cstarmap2Dlg

Cstarmap2Dlg::Cstarmap2Dlg(CWnd* pParent /*=nullptr*/)
    : CDialogEx(IDD_STARMAP2_DIALOG, pParent)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void Cstarmap2Dlg::DoDataExchange(CDataExchange* pDX)
{
    CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(Cstarmap2Dlg, CDialogEx)
    ON_BN_CLICKED(
        IDC_BUTTON_DISPLAY,
        &Cstarmap2Dlg::OnBnClickedButtonDisplay)
END_MESSAGE_MAP()


BOOL Cstarmap2Dlg::OnInitDialog()
{
    CDialogEx::OnInitDialog();

    SetIcon(m_hIcon, TRUE);
    SetIcon(m_hIcon, FALSE);

    SetDlgItemText(IDC_EDIT_ALT, _T("45.0"));
    SetDlgItemText(IDC_EDIT_AZ, _T("180.0"));

    return TRUE;
}


bool Cstarmap2Dlg::ReadCoordinate(
    int controlId,
    LPCTSTR label,
    double minimum,
    double maximum,
    bool includeMaximum,
    double& value) const
{
    CString text;
    GetDlgItemText(controlId, text);

    if (!ParseDoubleStrict(text, value))
    {
        CString message;
        message.Format(
            _T("%s must be a number."),
            label);

        AfxMessageBox(message);
        return false;
    }

    const bool inRange =
        includeMaximum
        ? (value >= minimum && value <= maximum)
        : (value >= minimum && value < maximum);

    if (!inRange)
    {
        CString message;

        if (includeMaximum)
        {
            message.Format(
                _T("%s must be between %.0f and %.0f degrees."),
                label,
                minimum,
                maximum);
        }
        else
        {
            message.Format(
                _T("%s must be between %.0f and less than %.0f degrees."),
                label,
                minimum,
                maximum);
        }

        AfxMessageBox(message);
        return false;
    }

    return true;
}


bool Cstarmap2Dlg::LaunchPythonStarMap(double altitude,double azimuth) const
{
    const CString starMapRoot =
        GetStarMapRoot();

    const CString pythonExe =
        GetPythonExePath();

    const CString mainPy =
        GetMainPyPath();

    if (GetFileAttributes(pythonExe) == INVALID_FILE_ATTRIBUTES)
    {
        AfxMessageBox(
            _T("Project Python environment was not found.\\n"
               "Please run Python\\setup_venv.ps1 first."));
        return false;
    }


    CString parameters;


    parameters.Format(
        _T("\"%s\" --altitude %.8f --azimuth %.8f"),
        mainPy.GetString(),
        altitude,
        azimuth
    );

    HINSTANCE result =
        ShellExecute(
            nullptr,
            _T("open"),
            pythonExe,
            parameters,
            starMapRoot,
            SW_SHOWNORMAL);

    if ((INT_PTR)result <= 32)
    {
        CString message;

        message.Format(
            _T("Unable to launch Python star map:\n%s"),
            pythonExe.GetString());

        AfxMessageBox(message);

        return false;
    }

    return true;
}


void Cstarmap2Dlg::OnBnClickedButtonDisplay()
{
    double altitude = 0.0;
    double azimuth = 0.0;

    if (!ReadCoordinate(
        IDC_EDIT_ALT,
        _T("Altitude"),
        -90.0,
        90.0,
        true,
        altitude))
    {
        return;
    }

    if (!ReadCoordinate(
        IDC_EDIT_AZ,
        _T("Azimuth"),
        0.0,
        360.0,
        false,
        azimuth))
    {
        return;
    }


    if (!LaunchPythonStarMap(altitude, azimuth))
    {
        return;
    }

    AfxMessageBox(
        _T("Telescope coordinates sent to Python FK6 star map."));
}
