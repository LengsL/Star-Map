// starmap2Dlg.cpp : implementation file
//

#include "pch.h"
#include "framework.h"
#include "starmap2.h"
#include "starmap2Dlg.h"
#include "afxdialogex.h"

#include <atlconv.h>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <Shellapi.h>
#include <tchar.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

namespace
{
    constexpr UINT_PTR GOTO_POLL_TIMER_ID = 1;

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

    CString CreateGotoResultPath()
    {
        TCHAR temporaryDirectory[MAX_PATH] = { 0 };
        if (GetTempPath(MAX_PATH, temporaryDirectory) == 0)
            return CString();

        CString path(temporaryDirectory);
        path.AppendFormat(_T("starmap2_goto_%lu_%llu.json"), GetCurrentProcessId(), GetTickCount64());
        return path;
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

    bool FindJsonNumber(const CString& document, LPCTSTR name, double& value)
    {
        CString key;
        key.Format(_T("\"%s\""), name);
        const int keyPosition = document.Find(key);
        if (keyPosition < 0)
            return false;

        const int valuePosition = document.Find(_T(':'), keyPosition + key.GetLength());
        if (valuePosition < 0)
            return false;

        const CString remainder = document.Mid(valuePosition + 1);
        const TCHAR* begin = remainder.GetString();
        TCHAR* end = nullptr;
        errno = 0;
        value = _tcstod(begin, &end);
        return end != begin && errno != ERANGE && std::isfinite(value);
    }

    bool FindJsonString(const CString& document, LPCTSTR name, CString& value)
    {
        CString key;
        key.Format(_T("\"%s\""), name);
        const int keyPosition = document.Find(key);
        if (keyPosition < 0)
            return false;
        const int colonPosition = document.Find(_T(':'), keyPosition + key.GetLength());
        const int openQuote = colonPosition < 0 ? -1 : document.Find(_T('"'), colonPosition + 1);
        if (openQuote < 0)
            return false;
        const int closeQuote = document.Find(_T('"'), openQuote + 1);
        if (closeQuote < 0)
            return false;
        value = document.Mid(openQuote + 1, closeQuote - openQuote - 1);
        return true;
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
    ON_WM_TIMER()
    ON_WM_DESTROY()
END_MESSAGE_MAP()


BOOL Cstarmap2Dlg::OnInitDialog()
{
    CDialogEx::OnInitDialog();

    SetIcon(m_hIcon, TRUE);
    SetIcon(m_hIcon, FALSE);

    SetDlgItemText(IDC_EDIT_ALT, _T("45.0"));
    SetDlgItemText(IDC_EDIT_AZ, _T("180.0"));
    SetDlgItemText(IDC_EDIT_LONGITUDE, _T("121.55"));
    SetDlgItemText(IDC_EDIT_LATITUDE, _T("29.87"));
    SetDlgItemText(IDC_EDIT_HEIGHT, _T("0.0"));
    SetDlgItemText(IDC_STATIC_GOTO_ALT, _T("Altitude: --"));
    SetDlgItemText(IDC_STATIC_GOTO_AZ, _T("Azimuth: --"));
    SetDlgItemText(IDC_STATIC_GOTO_STATUS, _T("No GOTO command received."));
    SetTimer(GOTO_POLL_TIMER_ID, 250, nullptr);

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
        const LPCTSTR unit = _tcscmp(label, _T("Height")) == 0 ? _T("metres") : _T("degrees");

        if (includeMaximum)
        {
            message.Format(
                _T("%s must be between %.0f and %.0f %s."),
                label,
                minimum,
                maximum,
                unit);
        }
        else
        {
            message.Format(
                _T("%s must be between %.0f and less than %.0f %s."),
                label,
                minimum,
                maximum,
                unit);
        }

        AfxMessageBox(message);
        return false;
    }

    return true;
}


bool Cstarmap2Dlg::LaunchPythonStarMap(
    double altitude,
    double azimuth,
    double longitude,
    double latitude,
    double height)
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

    m_gotoResultPath = CreateGotoResultPath();
    if (m_gotoResultPath.IsEmpty())
    {
        AfxMessageBox(_T("Unable to create a temporary GOTO result path."));
        return false;
    }
    DeleteFile(m_gotoResultPath);
    m_lastGotoSequence = 0;
    SetDlgItemText(IDC_STATIC_GOTO_ALT, _T("Altitude: --"));
    SetDlgItemText(IDC_STATIC_GOTO_AZ, _T("Azimuth: --"));
    SetDlgItemText(IDC_STATIC_GOTO_STATUS, _T("Waiting for a GOTO command."));

    parameters.Format(
        _T("\"%s\" --altitude %.8f --azimuth %.8f --longitude %.8f --latitude %.8f --height %.3f --goto-output \"%s\""),
        mainPy.GetString(),
        altitude,
        azimuth,
        longitude,
        latitude,
        height,
        m_gotoResultPath.GetString()
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
    double longitude = 0.0;
    double latitude = 0.0;
    double height = 0.0;

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

    if (!ReadCoordinate(
        IDC_EDIT_LONGITUDE,
        _T("Longitude"),
        -180.0,
        180.0,
        true,
        longitude) ||
        !ReadCoordinate(
            IDC_EDIT_LATITUDE,
            _T("Latitude"),
            -90.0,
            90.0,
            true,
            latitude) ||
        !ReadCoordinate(
            IDC_EDIT_HEIGHT,
            _T("Height"),
            -1000.0,
            100000.0,
            true,
            height))
    {
        return;
    }

    if (!LaunchPythonStarMap(altitude, azimuth, longitude, latitude, height))
    {
        return;
    }

    AfxMessageBox(
        _T("Actual telescope and observatory coordinates sent to the Python sky map."));
}


void Cstarmap2Dlg::ReadGotoResult()
{
    if (m_gotoResultPath.IsEmpty())
        return;

    CStdioFile file;
    if (!file.Open(m_gotoResultPath, CFile::modeRead | CFile::typeText | CFile::shareDenyNone))
        return;

    CString document;
    CString line;
    while (file.ReadString(line))
        document += line;

    double sequence = 0.0;
    double altitude = 0.0;
    double azimuth = 0.0;
    if (!FindJsonNumber(document, _T("sequence"), sequence) ||
        !FindJsonNumber(document, _T("goto_altitude_deg"), altitude) ||
        !FindJsonNumber(document, _T("goto_azimuth_deg"), azimuth) ||
        sequence < 0.0 || sequence != std::floor(sequence) ||
        sequence <= static_cast<double>(m_lastGotoSequence) ||
        altitude < -90.0 || altitude > 90.0 || azimuth < 0.0 || azimuth >= 360.0)
    {
        return;
    }

    m_lastGotoSequence = static_cast<ULONGLONG>(sequence);
    CString message;
    message.Format(_T("Altitude: %+.4f deg"), altitude);
    SetDlgItemText(IDC_STATIC_GOTO_ALT, message);
    message.Format(_T("Azimuth: %.4f deg"), azimuth);
    SetDlgItemText(IDC_STATIC_GOTO_AZ, message);

    CString objectName;
    CString command;
    FindJsonString(document, _T("object_name"), objectName);
    FindJsonString(document, _T("command"), command);
    message.Format(_T("%s command: %s"), command.GetString(), objectName.GetString());
    SetDlgItemText(IDC_STATIC_GOTO_STATUS, message);
}


void Cstarmap2Dlg::OnTimer(UINT_PTR nIDEvent)
{
    if (nIDEvent == GOTO_POLL_TIMER_ID)
        ReadGotoResult();
    CDialogEx::OnTimer(nIDEvent);
}


void Cstarmap2Dlg::OnDestroy()
{
    KillTimer(GOTO_POLL_TIMER_ID);
    CDialogEx::OnDestroy();
}
