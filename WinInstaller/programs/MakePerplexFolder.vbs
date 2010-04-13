' Make Working Copy.vbs - A script to copy the commonly-used 
'   Perplex programs and data files to a new folder on the
'   Desktop.
Option Explicit

Dim OpenAfterMaking, RunBuildAfterMaking

OpenAfterMaking = True
RunBuildAfterMaking = False

Dim DestinationFolder
DestinationFolder = BrowseFolder( "", False )

Dim Fso, MyFolder
Set Fso = CreateObject("Scripting.FileSystemObject")
MyFolder = fso.GetParentFolderName(wscript.ScriptFullName)

Call CopyIfNeeded("OpenCommandWindowHere.vbs", MyFolder, DestinationFolder)
Call CopyIfNeeded("What to do next.txt", MyFolder, DestinationFolder)
Call CopyIfNeeded("perplex_option.dat", MyFolder, DestinationFolder)
Call CopyIfNeeded("perplex_plot_option.dat", MyFolder, DestinationFolder)
Call CopyIfNeeded("hp02ver.dat", MyFolder+"\datafiles", DestinationFolder)
Call CopyIfNeeded("solution_models.dat", MyFolder+"\datafiles", DestinationFolder)

If OpenAfterMaking Then
	Dim Shell
	Set Shell = CreateObject("Shell.Application")
	Shell.Explore DestinationFolder
End If

If RunBuildAfterMaking Then
	Dim objShell
	Set objShell = CreateObject("WScript.Shell")
	objShell.Run "%comspec% /k cd "+DestinationFolder+ " & build.exe"
End If

Sub CopyIfNeeded(theItem, sourceFolder, destinationFolder)
' This function will copy a file (theItem) from a sourceFolder
' to a destinationFolder only if it exists and is not already
' present at the destination

'	WScript.Echo "Copying "+theItem+ " from " + sourceFolder + " to " + destinationFolder + "."
	Dim filesys
	Set filesys = CreateObject("Scripting.FileSystemObject")
	If filesys.FileExists(sourceFolder+"\"+theItem) Then
		If Not filesys.FileExists(destinationFolder+"\"+theItem) Then
			filesys.CopyFile sourceFolder+"\"+theItem, destinationFolder+"\"
		End If
	End If
End Sub


Function BrowseFolder( myStartLocation, blnSimpleDialog )
' This function generates a Browse Folder dialog
' and returns the selected folder as a string.
'
' Arguments:
' myStartLocation   [string]  start folder for dialog, or "My Computer", or
'                             empty string to open in "Desktop\My Documents"
' blnSimpleDialog   [boolean] if False, an additional text field will be
'                             displayed where the folder can be selected
'                             by typing the fully qualified path
'
' Returns:          [string]  the fully qualified path to the selected folder
'
' Based on the Hey Scripting Guys article
' "How Can I Show Users a Dialog Box That Only Lets Them Select Folders?"
' http://www.microsoft.com/technet/scriptcenter/resources/qanda/jun05/hey0617.mspx
'
' Function written by Rob van der Woude
' http://www.robvanderwoude.com
    Const MY_COMPUTER   = &H11&
    Const WINDOW_HANDLE = 0 ' Must ALWAYS be 0

    Dim numOptions, objFolder, objFolderItem
    Dim objPath, objShell, strPath, strPrompt

    ' Set the options for the dialog window
    strPrompt = "Select/create a folder to turn into a working copy of Perplex:"
    If blnSimpleDialog = True Then
        numOptions = 0      ' Simple dialog
    Else
        numOptions = &H10&  ' Additional text field to type folder path
    End If
    
    ' Create a Windows Shell object
    Set objShell = CreateObject( "Shell.Application" )

    ' If specified, convert "My Computer" to a valid
    ' path for the Windows Shell's BrowseFolder method
    If UCase( myStartLocation ) = "MY COMPUTER" Then
        Set objFolder = objShell.Namespace( MY_COMPUTER )
        Set objFolderItem = objFolder.Self
        strPath = objFolderItem.Path
    Else
        strPath = myStartLocation
    End If

    Set objFolder = objShell.BrowseForFolder( WINDOW_HANDLE, strPrompt, _
                                              numOptions, strPath )

    ' Quit if no folder was selected
    If objFolder Is Nothing Then
        BrowseFolder = ""
        Exit Function
    End If

    ' Retrieve the path of the selected folder
    Set objFolderItem = objFolder.Self
    objPath = objFolderItem.Path

    ' Return the path of the selected folder
    BrowseFolder = objPath
End Function