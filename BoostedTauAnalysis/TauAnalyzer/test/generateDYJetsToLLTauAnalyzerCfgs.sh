#!/bin/bash

if [ $# -gt 2 ]
    then
    echo "Usage: ./generateDYJetsToLLTauAnalyzerCfgs.sh <version> <template cfg>"
    exit 0
fi

####STUFF TO CONFIGURE####

#version
version=$1
templateCfg=$2
infoTag=""
dir=$version

#number of samples
nSamples=2
iBeg=0
iEnd=`expr $nSamples - 1`

#input file prefix
begInputFilePrefix="root://eoscms//eos/cms/store/user/friccita/DYJetsToLL_M-"
endInputFilePrefix="Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_v3/data_no_selection_"
inputFilePrefixM10To50="${begInputFilePrefix}10To50_TuneZ2Star_8TeV-madgraph-${endInputFilePrefix}"
inputFilePrefixM50="${begInputFilePrefix}50_TuneZ2Star_8TeV-madgraph-tarball-${endInputFilePrefix}"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/DYJetsToLL/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/DYJetsToLL/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefixM10To50}100_1_yns.root',\n    '${inputFilePrefixM10To50}101_1_Onx.root',\n    '${inputFilePrefixM10To50}102_1_Yf4.root',\n    '${inputFilePrefixM10To50}103_1_wwR.root',\n    '${inputFilePrefixM10To50}104_1_1Gy.root',\n    '${inputFilePrefixM10To50}105_1_wS4.root',\n    '${inputFilePrefixM10To50}106_1_6Ru.root',\n    '${inputFilePrefixM10To50}107_1_GP5.root',\n    '${inputFilePrefixM10To50}108_1_nWF.root',\n    '${inputFilePrefixM10To50}109_1_t0a.root',\n    '${inputFilePrefixM10To50}10_1_kmE.root',\n    '${inputFilePrefixM10To50}110_1_CVK.root',\n    '${inputFilePrefixM10To50}111_1_RuH.root',\n    '${inputFilePrefixM10To50}112_1_TKV.root',\n    '${inputFilePrefixM10To50}113_1_t27.root',\n    '${inputFilePrefixM10To50}114_1_1J1.root',\n    '${inputFilePrefixM10To50}115_1_ST0.root',\n    '${inputFilePrefixM10To50}116_1_QvR.root',\n    '${inputFilePrefixM10To50}117_1_fJu.root',\n    '${inputFilePrefixM10To50}118_1_IZt.root',\n    '${inputFilePrefixM10To50}119_1_WTc.root',\n    '${inputFilePrefixM10To50}11_1_rxB.root',\n    '${inputFilePrefixM10To50}120_1_Z8A.root',\n    '${inputFilePrefixM10To50}121_1_wgj.root',\n    '${inputFilePrefixM10To50}122_1_Nhy.root',\n    '${inputFilePrefixM10To50}123_1_bkH.root',\n    '${inputFilePrefixM10To50}124_1_0PT.root',\n    '${inputFilePrefixM10To50}125_1_IRm.root',\n    '${inputFilePrefixM10To50}126_1_2Ay.root',\n    '${inputFilePrefixM10To50}127_1_GFl.root',\n    '${inputFilePrefixM10To50}128_1_YGv.root',\n    '${inputFilePrefixM10To50}129_1_KYP.root',\n    '${inputFilePrefixM10To50}12_1_L8R.root',\n    '${inputFilePrefixM10To50}130_1_FQr.root',\n    '${inputFilePrefixM10To50}131_1_nPx.root',\n    '${inputFilePrefixM10To50}132_1_zON.root',\n    '${inputFilePrefixM10To50}133_1_J66.root',\n    '${inputFilePrefixM10To50}134_1_jvw.root',\n    '${inputFilePrefixM10To50}135_1_DRM.root',\n    '${inputFilePrefixM10To50}136_1_tqT.root',\n    '${inputFilePrefixM10To50}137_1_tkB.root',\n    '${inputFilePrefixM10To50}138_1_wmY.root',\n    '${inputFilePrefixM10To50}139_1_rf3.root',\n    '${inputFilePrefixM10To50}13_1_0Mc.root',\n    '${inputFilePrefixM10To50}140_1_pN5.root',\n    '${inputFilePrefixM10To50}141_1_Kvi.root',\n    '${inputFilePrefixM10To50}142_1_FD0.root',\n    '${inputFilePrefixM10To50}143_1_nTV.root',\n    '${inputFilePrefixM10To50}144_1_OqB.root',\n    '${inputFilePrefixM10To50}145_1_J4S.root',\n    '${inputFilePrefixM10To50}146_1_eC4.root',\n    '${inputFilePrefixM10To50}147_1_igd.root',\n    '${inputFilePrefixM10To50}148_1_pTt.root',\n    '${inputFilePrefixM10To50}149_1_FRK.root',\n    '${inputFilePrefixM10To50}14_1_dlM.root',\n    '${inputFilePrefixM10To50}150_1_0Os.root',\n    '${inputFilePrefixM10To50}151_1_CHQ.root',\n    '${inputFilePrefixM10To50}152_1_3aD.root',\n    '${inputFilePrefixM10To50}153_1_A8b.root',\n    '${inputFilePrefixM10To50}154_1_pap.root',\n    '${inputFilePrefixM10To50}155_1_nxn.root',\n    '${inputFilePrefixM10To50}156_1_qBb.root',\n    '${inputFilePrefixM10To50}157_1_zJ0.root',\n    '${inputFilePrefixM10To50}158_2_qXU.root',\n    '${inputFilePrefixM10To50}159_1_2Th.root',\n    '${inputFilePrefixM10To50}15_1_4N1.root',\n    '${inputFilePrefixM10To50}160_1_Ha6.root',\n    '${inputFilePrefixM10To50}161_1_RTY.root',\n    '${inputFilePrefixM10To50}162_1_Cdr.root',\n    '${inputFilePrefixM10To50}163_1_Rb1.root',\n    '${inputFilePrefixM10To50}16_1_8aL.root',\n    '${inputFilePrefixM10To50}17_1_qeJ.root',\n    '${inputFilePrefixM10To50}18_1_hbp.root',\n    '${inputFilePrefixM10To50}19_1_qYj.root',\n    '${inputFilePrefixM10To50}1_1_dDh.root',\n    '${inputFilePrefixM10To50}20_1_2qM.root',\n    '${inputFilePrefixM10To50}21_1_Ahb.root',\n    '${inputFilePrefixM10To50}22_1_sdi.root',\n    '${inputFilePrefixM10To50}23_1_ybd.root',\n    '${inputFilePrefixM10To50}24_1_CaB.root',\n    '${inputFilePrefixM10To50}25_1_ujp.root',\n    '${inputFilePrefixM10To50}26_1_4y5.root',\n    '${inputFilePrefixM10To50}27_1_6sh.root',\n    '${inputFilePrefixM10To50}28_1_TQN.root',\n    '${inputFilePrefixM10To50}29_1_9Ie.root',\n    '${inputFilePrefixM10To50}2_1_D0w.root',\n    '${inputFilePrefixM10To50}30_1_2Hr.root',\n    '${inputFilePrefixM10To50}31_1_4iY.root',\n    '${inputFilePrefixM10To50}32_1_v3e.root',\n    '${inputFilePrefixM10To50}33_1_tW3.root',\n    '${inputFilePrefixM10To50}34_1_Ea0.root',\n    '${inputFilePrefixM10To50}35_1_SzF.root',\n    '${inputFilePrefixM10To50}36_1_bO6.root',\n    '${inputFilePrefixM10To50}37_1_6Yc.root',\n    '${inputFilePrefixM10To50}38_1_X6x.root',\n    '${inputFilePrefixM10To50}39_1_kEO.root',\n    '${inputFilePrefixM10To50}3_1_Ks6.root',\n    '${inputFilePrefixM10To50}40_1_1H1.root',\n    '${inputFilePrefixM10To50}41_1_C4m.root',\n    '${inputFilePrefixM10To50}42_1_gx0.root',\n    '${inputFilePrefixM10To50}43_1_uTQ.root',\n    '${inputFilePrefixM10To50}44_1_Gu9.root',\n    '${inputFilePrefixM10To50}45_1_XwN.root',\n    '${inputFilePrefixM10To50}46_1_xbL.root',\n    '${inputFilePrefixM10To50}47_1_J3n.root',\n    '${inputFilePrefixM10To50}48_1_B8i.root',\n    '${inputFilePrefixM10To50}49_1_zSX.root',\n    '${inputFilePrefixM10To50}4_1_HtU.root',\n    '${inputFilePrefixM10To50}50_1_wLb.root',\n    '${inputFilePrefixM10To50}51_1_8Xh.root',\n    '${inputFilePrefixM10To50}52_1_yDL.root',\n    '${inputFilePrefixM10To50}53_1_7I3.root',\n    '${inputFilePrefixM10To50}54_1_kEG.root',\n    '${inputFilePrefixM10To50}55_1_3aW.root',\n    '${inputFilePrefixM10To50}56_1_z7i.root',\n    '${inputFilePrefixM10To50}57_1_eaf.root',\n    '${inputFilePrefixM10To50}58_1_hFS.root',\n    '${inputFilePrefixM10To50}59_1_trr.root',\n    '${inputFilePrefixM10To50}5_1_Ogw.root',\n    '${inputFilePrefixM10To50}60_1_qUt.root',\n    '${inputFilePrefixM10To50}61_1_t3n.root',\n    '${inputFilePrefixM10To50}62_1_qgz.root',\n    '${inputFilePrefixM10To50}63_1_ZCD.root',\n    '${inputFilePrefixM10To50}64_1_962.root',\n    '${inputFilePrefixM10To50}65_1_KBu.root',\n    '${inputFilePrefixM10To50}66_1_lbj.root',\n    '${inputFilePrefixM10To50}67_1_zzP.root',\n    '${inputFilePrefixM10To50}68_1_qMK.root',\n    '${inputFilePrefixM10To50}69_1_NyH.root',\n    '${inputFilePrefixM10To50}6_1_GOp.root',\n    '${inputFilePrefixM10To50}70_1_NwW.root',\n    '${inputFilePrefixM10To50}71_1_EGr.root',\n    '${inputFilePrefixM10To50}72_1_lKM.root',\n    '${inputFilePrefixM10To50}73_1_MdV.root',\n    '${inputFilePrefixM10To50}74_1_dyX.root',\n    '${inputFilePrefixM10To50}75_1_MG0.root',\n    '${inputFilePrefixM10To50}76_1_QMi.root',\n    '${inputFilePrefixM10To50}77_1_qQ1.root',\n    '${inputFilePrefixM10To50}78_1_yjA.root',\n    '${inputFilePrefixM10To50}79_1_87Q.root',\n    '${inputFilePrefixM10To50}7_1_SsT.root',\n    '${inputFilePrefixM10To50}80_1_OgL.root',\n    '${inputFilePrefixM10To50}81_1_8nF.root',\n    '${inputFilePrefixM10To50}82_1_r2B.root',\n    '${inputFilePrefixM10To50}83_1_QOa.root',\n    '${inputFilePrefixM10To50}84_1_JHx.root',\n    '${inputFilePrefixM10To50}85_1_1Xy.root',\n    '${inputFilePrefixM10To50}86_1_8Zk.root',\n    '${inputFilePrefixM10To50}87_1_HC5.root',\n    '${inputFilePrefixM10To50}88_1_yJ1.root',\n    '${inputFilePrefixM10To50}89_1_NXD.root',\n    '${inputFilePrefixM10To50}8_1_Orf.root',\n    '${inputFilePrefixM10To50}90_1_A0h.root',\n    '${inputFilePrefixM10To50}91_1_B0l.root',\n    '${inputFilePrefixM10To50}92_1_ljI.root',\n    '${inputFilePrefixM10To50}93_1_KJ0.root',\n    '${inputFilePrefixM10To50}94_1_X61.root',\n    '${inputFilePrefixM10To50}95_1_SCz.root',\n    '${inputFilePrefixM10To50}96_1_KEq.root',\n    '${inputFilePrefixM10To50}97_1_TqU.root',\n    '${inputFilePrefixM10To50}98_1_uGT.root',\n    '${inputFilePrefixM10To50}99_1_FyO.root',\n    '${inputFilePrefixM10To50}9_1_mFA.root'    ])" "readFiles.extend([\n    '${inputFilePrefixM50}100_1_Csf.root',\n    '${inputFilePrefixM50}101_1_J24.root',\n    '${inputFilePrefixM50}102_1_BUw.root',\n    '${inputFilePrefixM50}103_1_z8y.root',\n    '${inputFilePrefixM50}104_1_hDx.root',\n    '${inputFilePrefixM50}105_1_4q2.root',\n    '${inputFilePrefixM50}106_1_usU.root',\n    '${inputFilePrefixM50}107_1_kwp.root',\n    '${inputFilePrefixM50}108_1_zin.root',\n    '${inputFilePrefixM50}109_1_xo8.root',\n    '${inputFilePrefixM50}10_1_v4l.root',\n    '${inputFilePrefixM50}110_1_7IT.root',\n    '${inputFilePrefixM50}111_1_fRI.root',\n    '${inputFilePrefixM50}112_1_T3P.root',\n    '${inputFilePrefixM50}113_1_qb2.root',\n    '${inputFilePrefixM50}114_1_wdy.root',\n    '${inputFilePrefixM50}115_1_5q6.root',\n    '${inputFilePrefixM50}116_1_zvA.root',\n    '${inputFilePrefixM50}117_1_eU3.root',\n    '${inputFilePrefixM50}118_1_QFr.root',\n    '${inputFilePrefixM50}119_1_uPo.root',\n    '${inputFilePrefixM50}11_1_b5g.root',\n    '${inputFilePrefixM50}120_1_n6u.root',\n    '${inputFilePrefixM50}121_1_Ppb.root',\n    '${inputFilePrefixM50}122_1_zEw.root',\n    '${inputFilePrefixM50}123_1_dJJ.root',\n    '${inputFilePrefixM50}124_1_c49.root',\n    '${inputFilePrefixM50}125_1_PUV.root',\n    '${inputFilePrefixM50}126_1_Uw0.root',\n    '${inputFilePrefixM50}127_1_EMW.root',\n    '${inputFilePrefixM50}128_1_Smw.root',\n    '${inputFilePrefixM50}129_1_V2T.root',\n    '${inputFilePrefixM50}12_1_k0o.root',\n    '${inputFilePrefixM50}130_1_Z1R.root',\n    '${inputFilePrefixM50}131_2_JNE.root',\n    '${inputFilePrefixM50}132_1_0Ox.root',\n    '${inputFilePrefixM50}133_1_Fvb.root',\n    '${inputFilePrefixM50}134_1_sgl.root',\n    '${inputFilePrefixM50}135_1_J54.root',\n    '${inputFilePrefixM50}136_2_RlS.root',\n    '${inputFilePrefixM50}137_1_Muj.root',\n    '${inputFilePrefixM50}138_2_1s3.root',\n    '${inputFilePrefixM50}139_2_Sp9.root',\n    '${inputFilePrefixM50}13_1_aIa.root',\n    '${inputFilePrefixM50}140_1_D7q.root',\n    '${inputFilePrefixM50}141_2_WRV.root',\n    '${inputFilePrefixM50}142_2_LMj.root',\n    '${inputFilePrefixM50}143_1_BMb.root',\n    '${inputFilePrefixM50}144_1_45r.root',\n    '${inputFilePrefixM50}145_2_DH3.root',\n    '${inputFilePrefixM50}146_1_o0h.root',\n    '${inputFilePrefixM50}147_1_nMh.root',\n    '${inputFilePrefixM50}148_2_BX3.root',\n    '${inputFilePrefixM50}149_1_U9C.root',\n    '${inputFilePrefixM50}14_2_GrH.root',\n    '${inputFilePrefixM50}150_1_DjG.root',\n    '${inputFilePrefixM50}151_2_QHo.root',\n    '${inputFilePrefixM50}152_2_1Le.root',\n    '${inputFilePrefixM50}153_2_jty.root',\n    '${inputFilePrefixM50}154_2_vhP.root',\n    '${inputFilePrefixM50}155_2_V7n.root',\n    '${inputFilePrefixM50}156_2_lIk.root',\n    '${inputFilePrefixM50}157_2_8r5.root',\n    '${inputFilePrefixM50}158_2_ILn.root',\n    '${inputFilePrefixM50}159_1_Rzf.root',\n    '${inputFilePrefixM50}15_1_bL7.root',\n    '${inputFilePrefixM50}160_2_gKo.root',\n    '${inputFilePrefixM50}161_1_5Na.root',\n    '${inputFilePrefixM50}162_2_1Hr.root',\n    '${inputFilePrefixM50}163_1_Gtq.root',\n    '${inputFilePrefixM50}164_2_HKT.root',\n    '${inputFilePrefixM50}165_2_nPK.root',\n    '${inputFilePrefixM50}166_2_i1c.root',\n    '${inputFilePrefixM50}167_1_jtB.root',\n    '${inputFilePrefixM50}168_1_Eyl.root',\n    '${inputFilePrefixM50}169_1_JlD.root',\n    '${inputFilePrefixM50}16_1_ZAQ.root',\n    '${inputFilePrefixM50}170_1_O1Y.root',\n    '${inputFilePrefixM50}171_1_Qgn.root',\n    '${inputFilePrefixM50}172_1_DaL.root',\n    '${inputFilePrefixM50}173_1_yXh.root',\n    '${inputFilePrefixM50}174_1_cwW.root',\n    '${inputFilePrefixM50}175_1_UHL.root',\n    '${inputFilePrefixM50}176_1_2F9.root',\n    '${inputFilePrefixM50}177_1_Vzy.root',\n    '${inputFilePrefixM50}178_1_w4V.root',\n    '${inputFilePrefixM50}179_1_YZo.root',\n    '${inputFilePrefixM50}17_1_s6z.root',\n    '${inputFilePrefixM50}180_1_HFx.root',\n    '${inputFilePrefixM50}181_1_mPk.root',\n    '${inputFilePrefixM50}182_1_89l.root',\n    '${inputFilePrefixM50}183_1_5EL.root',\n    '${inputFilePrefixM50}184_1_FkX.root',\n    '${inputFilePrefixM50}185_1_NLz.root',\n    '${inputFilePrefixM50}186_1_gPD.root',\n    '${inputFilePrefixM50}187_1_Zto.root',\n    '${inputFilePrefixM50}188_1_Fuf.root',\n    '${inputFilePrefixM50}189_1_6Bv.root',\n    '${inputFilePrefixM50}18_1_lOH.root',\n    '${inputFilePrefixM50}190_1_QdJ.root',\n    '${inputFilePrefixM50}191_2_vfV.root',\n    '${inputFilePrefixM50}192_2_yq4.root',\n    '${inputFilePrefixM50}193_2_LMm.root',\n    '${inputFilePrefixM50}194_1_GdT.root',\n    '${inputFilePrefixM50}195_2_AKl.root',\n    '${inputFilePrefixM50}196_2_6ZP.root',\n    '${inputFilePrefixM50}197_1_tNT.root',\n    '${inputFilePrefixM50}198_1_j0G.root',\n    '${inputFilePrefixM50}199_2_MyX.root',\n    '${inputFilePrefixM50}19_1_eco.root',\n    '${inputFilePrefixM50}1_1_qTf.root',\n    '${inputFilePrefixM50}200_2_FS2.root',\n    '${inputFilePrefixM50}201_2_VPD.root',\n    '${inputFilePrefixM50}202_2_Ll7.root',\n    '${inputFilePrefixM50}203_2_quL.root',\n    '${inputFilePrefixM50}204_2_r4b.root',\n    '${inputFilePrefixM50}205_2_JsN.root',\n    '${inputFilePrefixM50}206_2_BOu.root',\n    '${inputFilePrefixM50}207_2_W04.root',\n    '${inputFilePrefixM50}208_2_tuW.root',\n    '${inputFilePrefixM50}209_1_7Wb.root',\n    '${inputFilePrefixM50}20_1_hIX.root',\n    '${inputFilePrefixM50}210_1_aaC.root',\n    '${inputFilePrefixM50}211_1_B2r.root',\n    '${inputFilePrefixM50}212_2_esT.root',\n    '${inputFilePrefixM50}213_2_VD8.root',\n    '${inputFilePrefixM50}214_1_RW2.root',\n    '${inputFilePrefixM50}215_1_PNj.root',\n    '${inputFilePrefixM50}216_1_5ma.root',\n    '${inputFilePrefixM50}217_1_N9Y.root',\n    '${inputFilePrefixM50}218_1_tGb.root',\n    '${inputFilePrefixM50}219_1_CyF.root',\n    '${inputFilePrefixM50}21_1_o3B.root',\n    '${inputFilePrefixM50}220_1_m7R.root',\n    '${inputFilePrefixM50}221_1_X3d.root',\n    '${inputFilePrefixM50}222_1_1X2.root',\n    '${inputFilePrefixM50}223_1_CZY.root',\n    '${inputFilePrefixM50}224_1_nL3.root',\n    '${inputFilePrefixM50}225_1_HtQ.root',\n    '${inputFilePrefixM50}226_1_ln8.root',\n    '${inputFilePrefixM50}227_1_LYp.root',\n    '${inputFilePrefixM50}228_1_nMa.root',\n    '${inputFilePrefixM50}229_1_RTs.root',\n    '${inputFilePrefixM50}22_2_KrR.root',\n    '${inputFilePrefixM50}230_1_8d0.root',\n    '${inputFilePrefixM50}231_1_RzL.root',\n    '${inputFilePrefixM50}232_1_mWT.root',\n    '${inputFilePrefixM50}233_1_faY.root',\n    '${inputFilePrefixM50}234_1_tpq.root',\n    '${inputFilePrefixM50}235_1_zGL.root',\n    '${inputFilePrefixM50}236_1_bOx.root',\n    '${inputFilePrefixM50}237_1_V8C.root',\n    '${inputFilePrefixM50}238_1_Blg.root',\n    '${inputFilePrefixM50}239_1_JCQ.root',\n    '${inputFilePrefixM50}23_1_m7Y.root',\n    '${inputFilePrefixM50}240_1_L2c.root',\n    '${inputFilePrefixM50}241_1_ul6.root',\n    '${inputFilePrefixM50}242_1_VOL.root',\n    '${inputFilePrefixM50}243_1_VjZ.root',\n    '${inputFilePrefixM50}244_1_qk3.root',\n    '${inputFilePrefixM50}245_1_kzE.root',\n    '${inputFilePrefixM50}246_1_zgh.root',\n    '${inputFilePrefixM50}247_1_9Sd.root',\n    '${inputFilePrefixM50}248_1_Qne.root',\n    '${inputFilePrefixM50}249_1_wCg.root',\n    '${inputFilePrefixM50}24_1_qXX.root',\n    '${inputFilePrefixM50}250_1_uC4.root',\n    '${inputFilePrefixM50}251_1_7rs.root',\n    '${inputFilePrefixM50}252_1_sOW.root',\n    '${inputFilePrefixM50}253_1_s5U.root',\n    '${inputFilePrefixM50}254_1_W2g.root',\n    '${inputFilePrefixM50}255_1_If9.root',\n    '${inputFilePrefixM50}256_1_Tfz.root',\n    '${inputFilePrefixM50}257_1_Rr7.root',\n    '${inputFilePrefixM50}258_1_9yb.root',\n    '${inputFilePrefixM50}259_1_znl.root',\n    '${inputFilePrefixM50}25_1_tTy.root',\n    '${inputFilePrefixM50}260_1_tFw.root',\n    '${inputFilePrefixM50}261_1_3aj.root',\n    '${inputFilePrefixM50}262_2_9Ux.root',\n    '${inputFilePrefixM50}263_1_gwf.root',\n    '${inputFilePrefixM50}264_1_94i.root',\n    '${inputFilePrefixM50}265_1_dRW.root',\n    '${inputFilePrefixM50}266_1_AdD.root',\n    '${inputFilePrefixM50}267_1_s6z.root',\n    '${inputFilePrefixM50}268_1_lRr.root',\n    '${inputFilePrefixM50}269_1_RvK.root',\n    '${inputFilePrefixM50}26_1_SAB.root',\n    '${inputFilePrefixM50}270_1_c7T.root',\n    '${inputFilePrefixM50}271_1_yJU.root',\n    '${inputFilePrefixM50}272_1_uhN.root',\n    '${inputFilePrefixM50}273_1_NCa.root',\n    '${inputFilePrefixM50}274_1_7T0.root',\n    '${inputFilePrefixM50}275_1_Ds1.root',\n    '${inputFilePrefixM50}276_1_VjS.root',\n    '${inputFilePrefixM50}277_1_7jT.root',\n    '${inputFilePrefixM50}278_1_hza.root',\n    '${inputFilePrefixM50}279_1_uKb.root',\n    '${inputFilePrefixM50}27_1_0ui.root',\n    '${inputFilePrefixM50}280_1_G3w.root',\n    '${inputFilePrefixM50}281_1_uY1.root',\n    '${inputFilePrefixM50}282_1_zD7.root',\n    '${inputFilePrefixM50}283_1_EM2.root',\n    '${inputFilePrefixM50}284_1_rMy.root',\n    '${inputFilePrefixM50}285_1_gh1.root',\n    '${inputFilePrefixM50}286_1_JeH.root',\n    '${inputFilePrefixM50}287_1_6QP.root',\n    '${inputFilePrefixM50}288_1_fUL.root',\n    '${inputFilePrefixM50}289_1_3oH.root',\n    '${inputFilePrefixM50}28_1_IHG.root',\n    '${inputFilePrefixM50}290_1_uCY.root',\n    '${inputFilePrefixM50}291_2_p6C.root',\n    '${inputFilePrefixM50}292_2_2IK.root',\n    '${inputFilePrefixM50}293_1_j6g.root',\n    '${inputFilePrefixM50}294_1_NNq.root',\n    '${inputFilePrefixM50}295_1_OW1.root'\n    ])\nreadFiles.extend([\n    '${inputFilePrefixM50}296_2_muk.root',\n    '${inputFilePrefixM50}297_2_IIX.root',\n    '${inputFilePrefixM50}298_2_rqS.root',\n    '${inputFilePrefixM50}299_1_KFg.root',\n    '${inputFilePrefixM50}29_1_ldx.root',\n    '${inputFilePrefixM50}2_1_wIP.root',\n    '${inputFilePrefixM50}300_2_8CX.root',\n    '${inputFilePrefixM50}301_2_iiV.root',\n    '${inputFilePrefixM50}302_2_3gl.root',\n    '${inputFilePrefixM50}303_2_ZFp.root',\n    '${inputFilePrefixM50}304_2_JUu.root',\n    '${inputFilePrefixM50}305_2_rqO.root',\n    '${inputFilePrefixM50}306_2_VKj.root',\n    '${inputFilePrefixM50}307_1_Jfi.root',\n    '${inputFilePrefixM50}308_2_psS.root',\n    '${inputFilePrefixM50}309_1_wMg.root',\n    '${inputFilePrefixM50}30_1_qHP.root',\n    '${inputFilePrefixM50}31_1_5AR.root',\n    '${inputFilePrefixM50}32_1_xHp.root',\n    '${inputFilePrefixM50}33_1_AvT.root',\n    '${inputFilePrefixM50}34_1_gYq.root',\n    '${inputFilePrefixM50}35_1_Fd1.root',\n    '${inputFilePrefixM50}36_1_QzZ.root',\n    '${inputFilePrefixM50}37_1_BoV.root',\n    '${inputFilePrefixM50}38_1_fX4.root',\n    '${inputFilePrefixM50}39_1_iMQ.root',\n    '${inputFilePrefixM50}3_1_WdR.root',\n    '${inputFilePrefixM50}40_1_Q7K.root',\n    '${inputFilePrefixM50}41_1_7vA.root',\n    '${inputFilePrefixM50}42_1_Fqf.root',\n    '${inputFilePrefixM50}43_1_mts.root',\n    '${inputFilePrefixM50}44_1_gMU.root',\n    '${inputFilePrefixM50}45_1_52f.root',\n    '${inputFilePrefixM50}46_1_nUc.root',\n    '${inputFilePrefixM50}47_1_c6s.root',\n    '${inputFilePrefixM50}48_1_WKY.root',\n    '${inputFilePrefixM50}49_1_p0t.root',\n    '${inputFilePrefixM50}4_1_lGZ.root',\n    '${inputFilePrefixM50}50_1_Un6.root',\n    '${inputFilePrefixM50}51_1_XgZ.root',\n    '${inputFilePrefixM50}52_1_bmi.root',\n    '${inputFilePrefixM50}53_1_1XP.root',\n    '${inputFilePrefixM50}54_1_RFv.root',\n    '${inputFilePrefixM50}55_1_M2c.root',\n    '${inputFilePrefixM50}56_1_evx.root',\n    '${inputFilePrefixM50}57_1_GFl.root',\n    '${inputFilePrefixM50}58_1_dbp.root',\n    '${inputFilePrefixM50}59_1_EPK.root',\n    '${inputFilePrefixM50}5_1_qg1.root',\n    '${inputFilePrefixM50}60_1_slF.root',\n    '${inputFilePrefixM50}61_1_2eQ.root',\n    '${inputFilePrefixM50}62_1_PXO.root',\n    '${inputFilePrefixM50}63_1_JhZ.root',\n    '${inputFilePrefixM50}64_1_gnW.root',\n    '${inputFilePrefixM50}65_1_EXB.root',\n    '${inputFilePrefixM50}66_1_uLO.root',\n    '${inputFilePrefixM50}67_1_2OY.root',\n    '${inputFilePrefixM50}68_1_9hN.root',\n    '${inputFilePrefixM50}69_1_vi6.root',\n    '${inputFilePrefixM50}6_1_Ete.root',\n    '${inputFilePrefixM50}70_1_Qax.root',\n    '${inputFilePrefixM50}71_1_p96.root',\n    '${inputFilePrefixM50}72_1_252.root',\n    '${inputFilePrefixM50}73_1_3mA.root',\n    '${inputFilePrefixM50}74_1_D3u.root',\n    '${inputFilePrefixM50}75_1_6Gp.root',\n    '${inputFilePrefixM50}76_1_2I3.root',\n    '${inputFilePrefixM50}77_1_BLm.root',\n    '${inputFilePrefixM50}78_1_eE9.root',\n    '${inputFilePrefixM50}79_1_XwD.root',\n    '${inputFilePrefixM50}7_1_DpY.root',\n    '${inputFilePrefixM50}80_1_XqT.root',\n    '${inputFilePrefixM50}81_1_Fbg.root',\n    '${inputFilePrefixM50}82_1_D6R.root',\n    '${inputFilePrefixM50}83_1_2NZ.root',\n    '${inputFilePrefixM50}84_1_yUb.root',\n    '${inputFilePrefixM50}85_1_IAZ.root',\n    '${inputFilePrefixM50}86_1_jUv.root',\n    '${inputFilePrefixM50}87_1_cD3.root',\n    '${inputFilePrefixM50}88_1_Yd7.root',\n    '${inputFilePrefixM50}89_1_ys8.root',\n    '${inputFilePrefixM50}8_1_XkG.root',\n    '${inputFilePrefixM50}90_1_Z0G.root',\n    '${inputFilePrefixM50}91_1_5Qs.root',\n    '${inputFilePrefixM50}92_1_BOO.root',\n    '${inputFilePrefixM50}93_1_TQw.root',\n    '${inputFilePrefixM50}94_1_KIo.root',\n    '${inputFilePrefixM50}95_1_anB.root',\n    '${inputFilePrefixM50}96_1_pMG.root',\n    '${inputFilePrefixM50}97_1_s7K.root',\n    '${inputFilePrefixM50}98_1_mri.root',\n    '${inputFilePrefixM50}99_1_7Tm.root',\n    '${inputFilePrefixM50}9_1_dS9.root'\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_DYJetsToLL_M-10To50.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_DYJetsToLL_M-50.root" )

#TauAnalyzer output files
highMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_DYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_DYJetsToLL_M-50_${version}.root" )
lowMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_DYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_DYJetsToLL_M-50_${version}.root" )
highMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_DYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_DYJetsToLL_M-50_${version}.root" )
lowMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_DYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_DYJetsToLL_M-50_${version}.root" )
highMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_DYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_DYJetsToLL_M-50_${version}.root" )
lowMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_DYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_DYJetsToLL_M-50_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}DYJetsToLL_M-10To50${infoTag}_${version}.root" "${EDMOutputFilePrefix}DYJetsToLL_M-50${infoTag}_${version}.root" )

#samples
samples=( "DYJetsToLL_M-10To50" "DYJetsToLL_M-50" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%HIGHMTNONISOTAUANALYZEROUTFILE%${highMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTALLTAUANALYZEROUTFILE%${highMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTISOTAUANALYZEROUTFILE%${highMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTNONISOTAUANALYZEROUTFILE%${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTALLTAUANALYZEROUTFILE%${lowMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTISOTAUANALYZEROUTFILE%${lowMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGGSREW%False%" -e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_cfg.py

  #generate job submission script for LSF
  cat <<EOF > tauanalyzer_${samples[${i}]}_cfg.sh
#!/bin/bash

jobDir="`pwd`"
fileNamePrefix="tauanalyzer_${samples[${i}]}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \$jobDir/\${fileNamePrefix}_cfg.py .
cmsRun \${fileNamePrefix}_cfg.py
cmsStage -f ${highMTNonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTNonIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${highMTIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTIsoTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${highMTAllTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
cmsStage -f ${lowMTAllTauAnalyzerOutputFiles[${i}]} /store/user/`whoami`/
rm ${highMTNonIsoTauAnalyzerOutputFiles[${i}]} ${lowMTNonIsoTauAnalyzerOutputFiles[${i}]} ${highMTIsoTauAnalyzerOutputFiles[${i}]} ${lowMTIsoTauAnalyzerOutputFiles[${i}]} ${highMTAllTauAnalyzerOutputFiles[${i}]} ${lowMTAllTauAnalyzerOutputFiles[${i}]}

exit 0
EOF
  chmod a+x tauanalyzer_${samples[${i}]}_cfg.sh
done

#generate run cfg that runs all files in the directory
cat <<EOF > runDYJetsToLLTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *DYJetsToLL*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runDYJetsToLLTauAnalyzerCfgs.sh

#generate script that submits all jobs to LSF
cat <<EOF > submitDYJetsToLLTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*DYJetsToLL*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 1nd -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitDYJetsToLLTauAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyDYJetsToLLFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "10To50" "50"
  do
  for cut in "Iso" "NonIso" ""
    do
    for MTBin in high low
      do
      cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_DYJetsToLL_M-\${sample}_${version}.root /data1/`whoami`/DYJetsToLL/analysis/
      cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_DYJetsToLL_M-\${sample}_${version}.root
    done
  done
done

exit 0
EOF
chmod a+x copyDYJetsToLLFromEOS.sh

exit 0
