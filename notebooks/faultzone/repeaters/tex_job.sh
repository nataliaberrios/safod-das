#!/bin/bash
#SBATCH -J tex
#SBATCH -p serc
#SBATCH -t 00:15:00
#SBATCH --mem=4G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/tex_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/tex_%j.err
ml system texlive 2>/dev/null || ml texlive 2>/dev/null
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
pdflatex -interaction=nonstopmode REPEATERS_progress.tex >/dev/null 2>&1
pdflatex -interaction=nonstopmode REPEATERS_progress.tex 2>&1 | tail -20
ls -l REPEATERS_progress.pdf
