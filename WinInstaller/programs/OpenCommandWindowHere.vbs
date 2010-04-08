' OpenCommandWindowHere.vbs - A script to start a new command window with the
'   working directory set to the script's location.
Option Explicit

Dim Fso, MyFolder
Set Fso = CreateObject("Scripting.FileSystemObject")
MyFolder = fso.GetParentFolderName(wscript.ScriptFullName)

	Dim objShell
	Set objShell = CreateObject("WScript.Shell")
	objShell.Run "%comspec% /k cd "+MyFolder+ " & echo To Get Started With Perplex, just type 'build'."

