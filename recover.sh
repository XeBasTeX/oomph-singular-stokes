#!/bin/bash
# Script pour sauvegarder 

rclone sync /home/bastien/Documents/oomph-lib-1.0.1331/user_drivers/NavierStokesInterns2019 remote:Sauvegarde/NavierStokesInterns2019 -v --filter-from .ignore
exit 0