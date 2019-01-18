#!/bin/bash

#   This script saves the results in temp-data to examples/data

#   NO ARGUMENTS

geom=$(awk 'BEGIN{FS=","; RS=""} {print $22}' temp-data/simulationParameters.csv)
nsites=$(awk 'BEGIN{FS=","; RS=""} {print $2}' temp-data/simulationParameters.csv)
dt_inv=$(awk 'BEGIN{FS=","; RS=""} {print $4}' temp-data/simulationParameters.csv)
beta=$(awk 'BEGIN{FS=","; RS=""} {print $6}' temp-data/simulationParameters.csv)
U=$(awk 'BEGIN{FS=","; RS=""} {print $12}' temp-data/simulationParameters.csv)
mu=$(awk 'BEGIN{FS=","; RS=""} {print $14}' temp-data/simulationParameters.csv)
ny=$(awk 'BEGIN{FS=","; RS=""} {print $24}' temp-data/simulationParameters.csv)

if [ $geom -eq 1 ]
then
  mkdir examples/data/1d-chain-pbc/
  mkdir examples/data/1d-chain-pbc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/1d-chain-pbc/
  mkdir examples/plots/1d-chain-pbc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/1d-chain-pbc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 2 ]
then
  mkdir examples/data/1d-chain-obc/
  mkdir examples/data/1d-chain-obc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/1d-chain-obc/
  mkdir examples/plots/1d-chain-obc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/1d-chain-obc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 3 ]
then
  mkdir examples/data/2d-sq-pbc/
  mkdir examples/data/2d-sq-pbc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/2d-sq-pbc/
  mkdir examples/plots/2d-sq-pbc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/2d-sq-pbc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 4 ]
then
  mkdir examples/data/2d-sq-obc/
  mkdir examples/data/2d-sq-obc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/2d-sq-obc/
  mkdir examples/plots/2d-sq-obc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/2d-sq-obc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 5 ]
then
  mkdir examples/data/2d-rec-pbc/
  mkdir examples/data/2d-rec-pbc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-rec-pbc/
  mkdir examples/plots/2d-rec-pbc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-rec-pbc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 6 ]
then
  mkdir examples/data/2d-rec-obc/
  mkdir examples/data/2d-rec-obc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-rec-obc/
  mkdir examples/plots/2d-rec-obc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-rec-obc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 7 ]
then
  mkdir examples/data/2d-triang-pbc/
  mkdir examples/data/2d-triang-pbc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-triang-pbc/
  mkdir examples/plots/2d-triang-pbc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-triang-pbc/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 8 ]
then
  mkdir examples/data/2d-triang-nanoribbon/
  mkdir examples/data/2d-triang-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-triang-nanoribbon/
  mkdir examples/plots/2d-triang-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-triang-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 9 ]
then
  mkdir examples/data/2d-hc-pbc/
  mkdir examples/data/2d-hc-pbc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/2d-hc-pbc/
  mkdir examples/plots/2d-hc-pbc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/2d-hc-pbc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 10 ]
then
  mkdir examples/data/2d-hc-nanoribbon/
  mkdir examples/data/2d-hc-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-hc-nanoribbon/
  mkdir examples/plots/2d-hc-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-hc-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 11 ]
then
  mkdir examples/data/2d-hc-nanoribbon-strain/
  mkdir examples/data/2d-hc-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-hc-nanoribbon-strain/
  mkdir examples/plots/2d-hc-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-hc-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 12 ]
then
  mkdir examples/data/2d-MoS2-pbc
  mkdir examples/data/2d-MoS2-pbc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/2d-MoS2-pbc
  mkdir examples/plots/2d-MoS2-pbc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/2d-MoS2-pbc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 13 ]
then
  mkdir examples/data/2d-MoS2-nanoribbon
  mkdir examples/data/2d-MoS2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-MoS2-nanoribbon
  mkdir examples/plots/2d-MoS2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-MoS2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 14 ]
then
  mkdir examples/data/2d-MoS2-nanoribbon-strain
  mkdir examples/data/2d-MoS2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-MoS2-nanoribbon-strain
  mkdir examples/plots/2d-MoS2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-MoS2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 15 ]
then
  mkdir examples/data/2d-WS2-pbc
  mkdir examples/data/2d-WS2-pbc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/2d-WS2-pbc
  mkdir examples/plots/2d-WS2-pbc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/2d-WS2-pbc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 16 ]
then
  mkdir examples/data/2d-WS2-nanoribbon
  mkdir examples/data/2d-WS2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-WS2-nanoribbon
  mkdir examples/plots/2d-WS2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-WS2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 17 ]
then
  mkdir examples/data/2d-WS2-nanoribbon-strain
  mkdir examples/data/2d-WS2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-WS2-nanoribbon-strain
  mkdir examples/plots/2d-WS2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-WS2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi


if [ $geom -eq 18 ]
then
  mkdir examples/data/2d-MoSe2-pbc
  mkdir examples/data/2d-MoSe2-pbc/N$nsites-BETA$beta-U$U-MU$mu
  mkdir examples/plots/2d-MoSe2-pbc
  mkdir examples/plots/2d-MoSe2-pbc/N$nsites-BETA$beta-U$U-MU$mu
  cp -r temp-data/* examples/data/2d-MoSe2-pbc/N$nsites-BETA$beta-U$U-MU$mu
fi

if [ $geom -eq 19 ]
then
  mkdir examples/data/2d-MoSe2-nanoribbon
  mkdir examples/data/2d-MoSe2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-MoSe2-nanoribbon
  mkdir examples/plots/2d-MoSe2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-MoSe2-nanoribbon/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi

if [ $geom -eq 20 ]
then
  mkdir examples/data/2d-MoSe2-nanoribbon-strain
  mkdir examples/data/2d-MoSe2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  mkdir examples/plots/2d-MoSe2-nanoribbon-strain
  mkdir examples/plots/2d-MoSe2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
  cp -r temp-data/* examples/data/2d-MoSe2-nanoribbon-strain/N$nsites-BETA$beta-U$U-MU$mu-NY$ny
fi
