
// starmap2.cpp : Defines the class behaviors for the application.
//

#include "pch.h"
#include "framework.h"
#include "starmap2.h"
#include "starmap2Dlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// Cstarmap2App

BEGIN_MESSAGE_MAP(Cstarmap2App, CWinApp)
END_MESSAGE_MAP()



Cstarmap2App::Cstarmap2App()
{
}


Cstarmap2App theApp;

BOOL Cstarmap2App::InitInstance()
{
	CWinApp::InitInstance();

	Cstarmap2Dlg dlg;
	m_pMainWnd = &dlg;

	dlg.DoModal();

	return FALSE;
}


