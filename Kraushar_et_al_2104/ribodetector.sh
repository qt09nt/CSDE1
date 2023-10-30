#!/bin/bash
#<------------------------Request for Resources----------------------->
#SBATCH -J ribodetector
#SBATCH --mem=30G
#SBATCH --time 01:00:00
#SBATCH -c 20
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1

