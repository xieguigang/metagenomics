
# copy files
/bin/cp -rf ./net5.0/* /usr/local/bin/
/bin/cp -rf ./16s.R /usr/local/bin/
/bin/cp -rf ./mothur_batch2.R /usr/local/bin/

# install package
R# --install.packages ./REnv.zip
R# --install.packages ./metagenomics.zip
R# --install.packages ./GCModeller-net5.0.zip