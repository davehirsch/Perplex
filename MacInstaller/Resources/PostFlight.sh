#!/bin/bash

# (Fails because the PATH in the installer-spawned process is not the PATH of the current user):
# Log my operation:
 echo "$SCRIPT_NAME: Preparing to test .bash_profile for the existence of location:$2"


# Create the .bash_profile file if necessary:
if [ ! -e $HOME/.bash_profile ];
then
	touch $HOME/.bash_profile
fi

# If the PATH doesn't currently contain the installed Perple_X files, then add that to the PATH
if ! grep -q "$2" $HOME/.bash_profile ;
then
	#append text to .bash_profile file - this will create the .bash_profile file if necessary
	echo "" >> $HOME/.bash_profile
	echo "###### Adding Perplex to default PATH " >> $HOME/.bash_profile
	echo "export PATH=\$PATH:$2" >> $HOME/.bash_profile
	echo "" >> $HOME/.bash_profile
fi

# If the PATH doesn't currently contain the installed Perple_X files, then add that to the PATH
if ! grep -q "\$DYLD_LIBRARY_PATH:/usr/local/lib/perplex" $HOME/.bash_profile ;
then
	#append text to .bash_profile file - this will create the .bash_profile file if necessary
	echo "" >> $HOME/.bash_profile
	echo "###### Adding  perplex libraries to DYLD path " >> $HOME/.bash_profile
	echo "export \$DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:/usr/local/lib/perplex" >> $HOME/.bash_profile
	echo "" >> $HOME/.bash_profile
fi

open -R "$2/MakePerplexFolder.app"
open "$2/MakePerplexFolder.app"
