
# copy files
/bin/cp -rf ./net5.0/* /usr/local/bin/
/bin/cp -rf ./16s.R /usr/local/bin/

# install package
R# --install.packages ./metagenomics.zip
R# --install.packages ./GCModeller_1.1.0-beta.zip