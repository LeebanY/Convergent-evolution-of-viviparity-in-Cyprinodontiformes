# EXAMPLE chaining script used for all 21 species.

#produce 2bit file.
faToTwoBit ../FASTA_files/GMassembly.fasta ../FASTA_files/GMassembly.fasta.2bit

# Produce chained file from the last alignment and the 2bit file.
axtChain -psl -linearGap=medium ../Goodeid_MAF/GM_XC.maf.psl ../FASTA_files/GMassembly.fasta.2bit ../FASTA_files/Xc_V1.0.fasta.2bit GM_XC.chain

# Chains with no chance of being netted are removed using pre-net chaining command. 
chainPreNet GM_XC.chain ../FASTA_files/GMassembly.fasta.sizes ../FASTA_files/Xc_V1.0.fasta.sizes GM_XC.prenet

# Chains are then netted and sorted.
chainNet GM_XC.prenet ../FASTA_files/GMassembly.fasta.sizes ../FASTA_files/Xc_V1.0.fasta.sizes GMXC.net XR.net
netSyntenic GMXC.net GMXC.syntenic.net
netToAxt GMXC.syntenic.net GM_XC.prenet ../FASTA_files/GMassembly.fasta.2bit ../FASTA_files/Xc_V1.0.fasta.2bit GM_netted_XC.axt
axtSort GM_netted_XC.axt GM_nettedsorted_XC.axt

# Chains are outputted then in maf format.
axtToMaf GM_nettedsorted_XC.axt ../FASTA_files/GMassembly.fasta.sizes ../FASTA_files/Xc_V1.0.fasta.sizes Girardinichthys_multiradiatus.Xenoophorus_captivus.sing.maf
