#!/bin/bash -l
#SBATCH -p regular    
#SBATCH -N 8        
#SBATCH -A m1542       
#SBATCH -t 00:15:00  
#SBATCH -J 2B2C    
/usr/bin/python main.py>log
