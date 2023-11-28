#! /bin/bash
sed -i "s/_/-/g" *.csv
sed -i "s/wB/ωB/g" *.csv
sed -i "s/wPBE/ωPBE/g" *.csv
