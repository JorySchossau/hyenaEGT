#!/bin/bash

rm -rf trial && ./explore.sh && ./collectend.sh trial && unzip -o trial.zip
