;----------------------------------------------------------------------------
;    MODULE NAME:   Perplex.MM
;
;        $Author:   USER "Dave Hirsch"  $
;      $Revision:   1.00  $
;          $Date:   26 Mar 2010 17:10:46  $
;
;----------------------------------------------------------------------------


;--- Include MAKEMSI support (with my customisations and MSI branding) ------
#define VER_FILENAME.VER  Perplex.VER      ;;I only want one VER file for all samples! (this line not actually required in "tryme.mm")
#include "Perplex.MMH"

;--- Define default location where file should install and add files --------

		<$DirectoryTree Key="DATADIR" Dir="c:\program files\Perple_X\datafiles" CHANGE="\" PrimaryFolder="N">
		<$DirectoryTree Key="CAIRODLLS" Dir="C:\WINDOWS\system32" CHANGE="\" PrimaryFolder="N">
		<$DirectoryTree Key="CODINGDLLS" Dir="C:\WINDOWS\system32" CHANGE="\" PrimaryFolder="N">
		<$DirectoryTree Key="INSTALLDIR" Dir="c:\program files\Perple_X" CHANGE="\" PrimaryFolder="Y">

<$Feature "MainFiles" Title=^Main Program Files^ Description="The main program files, such as vertex.exe, which are the core of Perplex.">
	<$Component "MainFiles" Create="Y" Directory_="INSTALLDIR">
		<$FilesExclude ".\programs\MakePerplexFolder.vbs" EXLIST="NotThese">   ;;Don't include these!
		<$FilesExclude ".\programs\Read Me.rtf" EXLIST="NotThese">   ;;Don't include these!
		<$FilesExclude ".\programs\Getting Started.rtf" EXLIST="NotThese">   ;;Don't include these!
		<$Files ".\programs\*" KeyFile="vertex.exe" EXLIST="NotThese">
		<$Path   "[INSTALLDIR]">
	<$/Component>

	<$Icon "perplex.ico" Key="@PerplexIcon">
	<$Component "ShortCutFile" Create="Y" Directory_="INSTALLDIR">
		;--- Create directory (add to "RemoveFile" table to keep validations happy) ---
		<$DirectoryTree Key="MY_SHORTCUT_FOLDER" Dir="[ProgramMenuFolder]\PerpleX" MAKE="Y" REMOVE=Y">

		;--- Add the file (make the component's keypath) ------------------------
		<$File RowKey="MWC.vbs" Source=".\programs\MakePerplexFolder.vbs" Comment='Simple tool to display environmental version information.' KeyPath="Y">
		<$File RowKey="ReadMe.rtf" Source=".\programs\Read Me.rtf" Comment='Basic Perplex Information'>
		<$File RowKey="Start.rtf" Source=".\programs\Getting Started.rtf" Comment='Basic Perplex Information'>

		;--- Create the shortcut (non-advertised) --------------------------------
	   #(
		   <$Shortcut
					   Dir="MY_SHORTCUT_FOLDER"
					Target="[!MWC.vbs]"
					 Title="MWC_Pplx|Make Working Copy of Perplex Folder"
			   Description="Makes a copy of the Perplex folder with commonly used programs and data files."
					  Icon="@PerplexIcon"
				   WorkDir="INSTALLDIR"
		   >
		   <$Shortcut
					   Dir="MY_SHORTCUT_FOLDER"
					Target="[!ReadMe.rtf]"
					 Title="Read Me"
			   Description="Basic Perplex Installation Information."
				   WorkDir="INSTALLDIR"
		   >
		   <$Shortcut
					   Dir="MY_SHORTCUT_FOLDER"
					Target="[!Start.rtf]"
					 Title="Getting Started with Perplex"
			   Description="Basic Perplex Tutorial."
				   WorkDir="INSTALLDIR"
		   >
	   #)

		;--- As mentioned above we can ignore validation warnings as install is per-machine ---
		<$MsiValFilter "Component ShortCutFile has non-advertised shortcuts. It should use a registry key under HKCU as its KeyPath"  Comment=^Not true, this MSI installs per-machine^>
		<$MsiValFilter "Component 'ShortCutFile' has both per-user and per-machine data with a per-machine KeyPath."                  Comment=^Not true, this MSI installs per-machine^>
	<$/Component>
<$/Feature>

<$Feature "DataFiles" Title=^Data Files^ Description="The set of data files required by Perplex (e.g., hp02ver.dat)">
	<$Component "DataFiles" Create="Y" Directory_="DATADIR">
		<$Files ".\datafiles\*">
	<$/Component>
<$/Feature>

<$Feature "Libraries" Title=^Graphics Libraries (DLLs)^ Description="These are libraries on which the graphics portions of Perple_X rely.  They include Cairo and its dependents.  You should either install Cairo yourself (http://cairographics.org) or install these libraries.">
	<$Component "CAIRODLLS" Create="Y" Directory_="CAIRODLLS">
		<$Files ".\graphics_dlls\*" KeyFile="libcairo-2.dll">
	<$/Component>
<$/Feature>

<$Feature "ReqLibraries" Title=^Required Libraries (DLLs)^ Description="These are programming libraries (C, Fortran) on which core portions of Perple_X rely. You should either install them yourself or install these libraries.">
	<$Component "CODINGDLLS" Create="Y" Directory_="CODINGDLLS">
		<$Files ".\programming_dlls\*" KeyFile="libgfortran-3.dll">
	<$/Component>
<$/Feature>

