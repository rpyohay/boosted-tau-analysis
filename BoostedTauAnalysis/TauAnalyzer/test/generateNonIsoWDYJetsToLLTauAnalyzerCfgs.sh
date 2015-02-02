#!/bin/bash

if [ $# -gt 2 ]
    then
    echo "Usage: ./generateNonIsoWDYJetsToLLTauAnalyzerCfgs.sh <version> <template cfg>"
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
endInputFilePrefix="Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_nonIsoW/data_no_selection"
inputFilePrefixM10To50="${begInputFilePrefix}10To50_TuneZ2Star_8TeV-madgraph-${endInputFilePrefix}"
inputFilePrefixM50="${begInputFilePrefix}50_TuneZ2Star_8TeV-madgraph-tarball-${endInputFilePrefix}"

#CleanJets output file prefix
#cleanJetsOutputFilePrefix="`pwd`/${dir}/"
cleanJetsOutputFilePrefix=""

#TauAnalyzer output file prefix
#tauAnalyzerOutputFilePrefix="/data1/yohay/DYJetsToLL/analysis/"
tauAnalyzerOutputFilePrefix=""

#EDM output file prefix
EDMOutputFilePrefix="/data1/`whoami`/nonIsoWDYJetsToLL/EDM_files/"

####VECTORS OF QUANTITIES FOR EACH SAMPLE####

#vector of input file blocks for each sample
#"readFiles.extend([\n    ])"
#    ])" "readFiles.extend([\n#
inputFileBlocks=( "readFiles.extend([\n    '${inputFilePrefixM10To50}_100_3_rRZ.root',\n    '${inputFilePrefixM10To50}_101_3_Z12.root',\n    '${inputFilePrefixM10To50}_102_2_SjT.root',\n    '${inputFilePrefixM10To50}_103_1_lO9.root',\n    '${inputFilePrefixM10To50}_104_1_XlH.root',\n    '${inputFilePrefixM10To50}_105_2_P8s.root',\n    '${inputFilePrefixM10To50}_106_1_7JY.root',\n    '${inputFilePrefixM10To50}_107_2_3a1.root',\n    '${inputFilePrefixM10To50}_108_1_22c.root',\n    '${inputFilePrefixM10To50}_109_3_ivD.root',\n    '${inputFilePrefixM10To50}_10_3_qLE.root',\n    '${inputFilePrefixM10To50}_110_3_wHa.root',\n    '${inputFilePrefixM10To50}_111_2_Ts8.root',\n    '${inputFilePrefixM10To50}_112_2_uB1.root',\n    '${inputFilePrefixM10To50}_113_3_QsY.root',\n    '${inputFilePrefixM10To50}_114_1_QDx.root',\n    '${inputFilePrefixM10To50}_115_1_uQ8.root',\n    '${inputFilePrefixM10To50}_116_1_Zr3.root',\n    '${inputFilePrefixM10To50}_117_3_pWX.root',\n    '${inputFilePrefixM10To50}_118_1_4wX.root',\n    '${inputFilePrefixM10To50}_119_1_jUa.root',\n    '${inputFilePrefixM10To50}_11_1_m9a.root',\n    '${inputFilePrefixM10To50}_120_1_ieL.root',\n    '${inputFilePrefixM10To50}_121_1_wDf.root',\n    '${inputFilePrefixM10To50}_122_1_Rwk.root',\n    '${inputFilePrefixM10To50}_123_1_n9X.root',\n    '${inputFilePrefixM10To50}_124_1_n0i.root',\n    '${inputFilePrefixM10To50}_125_1_aPe.root',\n    '${inputFilePrefixM10To50}_126_1_MqR.root',\n    '${inputFilePrefixM10To50}_127_1_pBN.root',\n    '${inputFilePrefixM10To50}_128_1_j7M.root',\n    '${inputFilePrefixM10To50}_129_1_GVU.root',\n    '${inputFilePrefixM10To50}_12_1_xvI.root',\n    '${inputFilePrefixM10To50}_130_1_Bk9.root',\n    '${inputFilePrefixM10To50}_131_1_Eb5.root',\n    '${inputFilePrefixM10To50}_132_1_2Nw.root',\n    '${inputFilePrefixM10To50}_133_1_ukH.root',\n    '${inputFilePrefixM10To50}_134_1_Irt.root',\n    '${inputFilePrefixM10To50}_135_1_QeI.root',\n    '${inputFilePrefixM10To50}_136_1_8A0.root',\n    '${inputFilePrefixM10To50}_137_1_ePa.root',\n    '${inputFilePrefixM10To50}_138_3_HMp.root',\n    '${inputFilePrefixM10To50}_139_1_4LY.root',\n    '${inputFilePrefixM10To50}_13_3_XFA.root',\n    '${inputFilePrefixM10To50}_140_1_g78.root',\n    '${inputFilePrefixM10To50}_141_3_Hpi.root',\n    '${inputFilePrefixM10To50}_142_3_JUc.root',\n    '${inputFilePrefixM10To50}_143_1_qbN.root',\n    '${inputFilePrefixM10To50}_144_1_3rD.root',\n    '${inputFilePrefixM10To50}_145_1_XuJ.root',\n    '${inputFilePrefixM10To50}_146_1_bUD.root',\n    '${inputFilePrefixM10To50}_147_1_t5Y.root',\n    '${inputFilePrefixM10To50}_148_1_ORz.root',\n    '${inputFilePrefixM10To50}_149_1_sNN.root',\n    '${inputFilePrefixM10To50}_14_3_eJt.root',\n    '${inputFilePrefixM10To50}_150_1_KP1.root',\n    '${inputFilePrefixM10To50}_151_1_R5F.root',\n    '${inputFilePrefixM10To50}_152_2_OoZ.root',\n    '${inputFilePrefixM10To50}_153_1_0Nl.root',\n    '${inputFilePrefixM10To50}_154_1_VJe.root',\n    '${inputFilePrefixM10To50}_155_3_E5I.root',\n    '${inputFilePrefixM10To50}_156_1_BFU.root',\n    '${inputFilePrefixM10To50}_157_2_pxc.root',\n    '${inputFilePrefixM10To50}_158_3_g1Y.root',\n    '${inputFilePrefixM10To50}_159_1_Qdg.root',\n    '${inputFilePrefixM10To50}_15_4_CMm.root',\n    '${inputFilePrefixM10To50}_160_1_rjA.root',\n    '${inputFilePrefixM10To50}_161_2_r5b.root',\n    '${inputFilePrefixM10To50}_162_1_mnH.root',\n    '${inputFilePrefixM10To50}_163_1_lrW.root',\n    '${inputFilePrefixM10To50}_16_1_spA.root',\n    '${inputFilePrefixM10To50}_17_3_Dxx.root',\n    '${inputFilePrefixM10To50}_18_4_iBe.root',\n    '${inputFilePrefixM10To50}_19_3_Hc3.root',\n    '${inputFilePrefixM10To50}_1_1_fJ3.root',\n    '${inputFilePrefixM10To50}_20_3_VSW.root',\n    '${inputFilePrefixM10To50}_21_2_DEM.root',\n    '${inputFilePrefixM10To50}_22_1_xry.root',\n    '${inputFilePrefixM10To50}_23_1_nas.root',\n    '${inputFilePrefixM10To50}_24_1_Yxp.root',\n    '${inputFilePrefixM10To50}_25_1_dku.root',\n    '${inputFilePrefixM10To50}_26_2_ifb.root',\n    '${inputFilePrefixM10To50}_27_1_Zkn.root',\n    '${inputFilePrefixM10To50}_28_1_ubh.root',\n    '${inputFilePrefixM10To50}_29_3_n3r.root',\n    '${inputFilePrefixM10To50}_2_2_6Jr.root',\n    '${inputFilePrefixM10To50}_30_4_kaa.root',\n    '${inputFilePrefixM10To50}_31_2_lv0.root',\n    '${inputFilePrefixM10To50}_32_1_mwh.root',\n    '${inputFilePrefixM10To50}_33_1_VJB.root',\n    '${inputFilePrefixM10To50}_34_2_W7d.root',\n    '${inputFilePrefixM10To50}_35_2_RsV.root',\n    '${inputFilePrefixM10To50}_36_2_4De.root',\n    '${inputFilePrefixM10To50}_37_3_MEJ.root',\n    '${inputFilePrefixM10To50}_38_2_pOK.root',\n    '${inputFilePrefixM10To50}_39_1_qUP.root',\n    '${inputFilePrefixM10To50}_3_1_c7B.root',\n    '${inputFilePrefixM10To50}_40_2_YpA.root',\n    '${inputFilePrefixM10To50}_41_3_M6u.root',\n    '${inputFilePrefixM10To50}_42_2_6pW.root',\n    '${inputFilePrefixM10To50}_43_3_uBp.root',\n    '${inputFilePrefixM10To50}_44_1_2Ll.root',\n    '${inputFilePrefixM10To50}_45_2_z5S.root',\n    '${inputFilePrefixM10To50}_46_3_muG.root',\n    '${inputFilePrefixM10To50}_47_3_iJk.root',\n    '${inputFilePrefixM10To50}_48_2_fxg.root',\n    '${inputFilePrefixM10To50}_49_2_RAR.root',\n    '${inputFilePrefixM10To50}_4_2_6w1.root',\n    '${inputFilePrefixM10To50}_50_2_0d6.root',\n    '${inputFilePrefixM10To50}_51_3_NJg.root',\n    '${inputFilePrefixM10To50}_52_2_dzn.root',\n    '${inputFilePrefixM10To50}_53_1_axh.root',\n    '${inputFilePrefixM10To50}_54_1_a65.root',\n    '${inputFilePrefixM10To50}_55_1_N9g.root',\n    '${inputFilePrefixM10To50}_56_1_DK4.root',\n    '${inputFilePrefixM10To50}_57_1_fwP.root',\n    '${inputFilePrefixM10To50}_58_3_efn.root',\n    '${inputFilePrefixM10To50}_59_1_14V.root',\n    '${inputFilePrefixM10To50}_5_3_9j8.root',\n    '${inputFilePrefixM10To50}_60_1_yHT.root',\n    '${inputFilePrefixM10To50}_61_2_RUg.root',\n    '${inputFilePrefixM10To50}_62_1_xT4.root',\n    '${inputFilePrefixM10To50}_63_1_d96.root',\n    '${inputFilePrefixM10To50}_64_3_exD.root',\n    '${inputFilePrefixM10To50}_65_3_yXo.root',\n    '${inputFilePrefixM10To50}_66_3_4w1.root',\n    '${inputFilePrefixM10To50}_67_1_01l.root',\n    '${inputFilePrefixM10To50}_68_2_8qu.root',\n    '${inputFilePrefixM10To50}_69_1_s0v.root',\n    '${inputFilePrefixM10To50}_6_3_mhs.root',\n    '${inputFilePrefixM10To50}_70_2_Y80.root',\n    '${inputFilePrefixM10To50}_71_1_Djc.root',\n    '${inputFilePrefixM10To50}_72_1_0OC.root',\n    '${inputFilePrefixM10To50}_73_3_ECV.root',\n    '${inputFilePrefixM10To50}_74_1_DSK.root',\n    '${inputFilePrefixM10To50}_75_2_MFY.root',\n    '${inputFilePrefixM10To50}_76_1_o6F.root',\n    '${inputFilePrefixM10To50}_77_3_GgD.root',\n    '${inputFilePrefixM10To50}_78_2_ite.root',\n    '${inputFilePrefixM10To50}_79_3_ckw.root',\n    '${inputFilePrefixM10To50}_7_3_haR.root',\n    '${inputFilePrefixM10To50}_80_1_i9P.root',\n    '${inputFilePrefixM10To50}_81_1_WKg.root',\n    '${inputFilePrefixM10To50}_82_3_61o.root',\n    '${inputFilePrefixM10To50}_83_3_7qq.root',\n    '${inputFilePrefixM10To50}_84_2_vPi.root',\n    '${inputFilePrefixM10To50}_85_1_yVO.root',\n    '${inputFilePrefixM10To50}_86_1_92b.root',\n    '${inputFilePrefixM10To50}_87_1_GTk.root',\n    '${inputFilePrefixM10To50}_88_3_cdG.root',\n    '${inputFilePrefixM10To50}_89_3_hJr.root',\n    '${inputFilePrefixM10To50}_8_3_Oro.root',\n    '${inputFilePrefixM10To50}_90_3_BVK.root',\n    '${inputFilePrefixM10To50}_91_1_oaf.root',\n    '${inputFilePrefixM10To50}_92_3_lC3.root',\n    '${inputFilePrefixM10To50}_93_2_Pwy.root',\n    '${inputFilePrefixM10To50}_94_1_7Ud.root',\n    '${inputFilePrefixM10To50}_95_1_gGh.root',\n    '${inputFilePrefixM10To50}_96_1_TOI.root',\n    '${inputFilePrefixM10To50}_97_1_EvR.root',\n    '${inputFilePrefixM10To50}_98_1_XpE.root',\n    '${inputFilePrefixM10To50}_99_1_KHv.root',\n    '${inputFilePrefixM10To50}_9_1_lgq.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefixM50}_100_1_n4a.root',\n    '${inputFilePrefixM50}_101_1_Fbx.root',\n    '${inputFilePrefixM50}_102_1_IEi.root',\n    '${inputFilePrefixM50}_103_1_ks2.root',\n    '${inputFilePrefixM50}_104_2_FBr.root',\n    '${inputFilePrefixM50}_105_1_Vkj.root',\n    '${inputFilePrefixM50}_106_1_9oz.root',\n    '${inputFilePrefixM50}_107_1_SlO.root',\n    '${inputFilePrefixM50}_108_1_VQX.root',\n    '${inputFilePrefixM50}_109_1_yGO.root',\n    '${inputFilePrefixM50}_10_1_GrB.root',\n    '${inputFilePrefixM50}_110_1_z01.root',\n    '${inputFilePrefixM50}_111_1_VBA.root',\n    '${inputFilePrefixM50}_112_1_PTv.root',\n    '${inputFilePrefixM50}_113_1_xV2.root',\n    '${inputFilePrefixM50}_114_1_JSd.root',\n    '${inputFilePrefixM50}_115_1_ZD1.root',\n    '${inputFilePrefixM50}_116_1_WQQ.root',\n    '${inputFilePrefixM50}_117_1_Hru.root',\n    '${inputFilePrefixM50}_118_1_MPj.root',\n    '${inputFilePrefixM50}_119_1_Ryc.root',\n    '${inputFilePrefixM50}_11_1_1zF.root',\n    '${inputFilePrefixM50}_120_1_9PZ.root',\n    '${inputFilePrefixM50}_121_1_o4t.root',\n    '${inputFilePrefixM50}_122_1_OGN.root',\n    '${inputFilePrefixM50}_123_1_AES.root',\n    '${inputFilePrefixM50}_124_1_P1O.root',\n    '${inputFilePrefixM50}_125_1_pDI.root',\n    '${inputFilePrefixM50}_126_1_Q5c.root',\n    '${inputFilePrefixM50}_127_1_XTS.root',\n    '${inputFilePrefixM50}_128_1_ZC6.root',\n    '${inputFilePrefixM50}_129_1_Wx3.root',\n    '${inputFilePrefixM50}_12_1_Wig.root',\n    '${inputFilePrefixM50}_130_1_SUf.root',\n    '${inputFilePrefixM50}_131_1_SCA.root',\n    '${inputFilePrefixM50}_132_1_UBd.root',\n    '${inputFilePrefixM50}_133_1_g8Z.root',\n    '${inputFilePrefixM50}_134_1_B4F.root',\n    '${inputFilePrefixM50}_135_1_47S.root',\n    '${inputFilePrefixM50}_136_1_qDd.root',\n    '${inputFilePrefixM50}_137_1_UHX.root',\n    '${inputFilePrefixM50}_138_1_qUW.root',\n    '${inputFilePrefixM50}_139_1_APV.root',\n    '${inputFilePrefixM50}_13_1_5vY.root',\n    '${inputFilePrefixM50}_140_1_9cg.root',\n    '${inputFilePrefixM50}_141_1_zI9.root',\n    '${inputFilePrefixM50}_142_1_tMI.root',\n    '${inputFilePrefixM50}_143_1_hrT.root',\n    '${inputFilePrefixM50}_144_1_5gE.root',\n    '${inputFilePrefixM50}_145_1_4Wz.root',\n    '${inputFilePrefixM50}_146_1_dhb.root',\n    '${inputFilePrefixM50}_147_1_6NL.root',\n    '${inputFilePrefixM50}_148_1_K8J.root',\n    '${inputFilePrefixM50}_149_1_Lxs.root',\n    '${inputFilePrefixM50}_14_1_YgI.root',\n    '${inputFilePrefixM50}_150_1_cIj.root',\n    '${inputFilePrefixM50}_151_1_lf4.root',\n    '${inputFilePrefixM50}_152_1_tTO.root',\n    '${inputFilePrefixM50}_153_1_vd9.root',\n    '${inputFilePrefixM50}_154_1_ecR.root',\n    '${inputFilePrefixM50}_155_1_Nux.root',\n    '${inputFilePrefixM50}_156_1_2ze.root',\n    '${inputFilePrefixM50}_157_1_YxJ.root',\n    '${inputFilePrefixM50}_158_1_Wfx.root',\n    '${inputFilePrefixM50}_159_1_StX.root',\n    '${inputFilePrefixM50}_15_1_cjg.root',\n    '${inputFilePrefixM50}_160_1_l6k.root',\n    '${inputFilePrefixM50}_161_1_Gwn.root',\n    '${inputFilePrefixM50}_162_1_rHe.root',\n    '${inputFilePrefixM50}_163_1_foD.root',\n    '${inputFilePrefixM50}_164_1_WRq.root',\n    '${inputFilePrefixM50}_165_1_mLy.root',\n    '${inputFilePrefixM50}_166_1_2bR.root',\n    '${inputFilePrefixM50}_167_1_OI2.root',\n    '${inputFilePrefixM50}_168_1_zoh.root',\n    '${inputFilePrefixM50}_169_1_AqG.root',\n    '${inputFilePrefixM50}_16_1_QDR.root',\n    '${inputFilePrefixM50}_170_1_jmn.root',\n    '${inputFilePrefixM50}_171_1_Xq8.root',\n    '${inputFilePrefixM50}_172_1_BnK.root',\n    '${inputFilePrefixM50}_173_1_ghK.root',\n    '${inputFilePrefixM50}_174_1_ZEn.root',\n    '${inputFilePrefixM50}_175_1_LN8.root',\n    '${inputFilePrefixM50}_176_1_Gqf.root',\n    '${inputFilePrefixM50}_177_1_ueq.root',\n    '${inputFilePrefixM50}_178_1_RBi.root',\n    '${inputFilePrefixM50}_179_1_5UE.root',\n    '${inputFilePrefixM50}_17_1_Qxy.root',\n    '${inputFilePrefixM50}_180_1_0q1.root',\n    '${inputFilePrefixM50}_181_1_TXc.root',\n    '${inputFilePrefixM50}_182_1_DWO.root',\n    '${inputFilePrefixM50}_183_1_j0i.root',\n    '${inputFilePrefixM50}_184_1_JRv.root',\n    '${inputFilePrefixM50}_185_1_Cdt.root',\n    '${inputFilePrefixM50}_186_1_3mP.root',\n    '${inputFilePrefixM50}_187_1_qUN.root',\n    '${inputFilePrefixM50}_188_1_AxZ.root',\n    '${inputFilePrefixM50}_189_1_Ms0.root',\n    '${inputFilePrefixM50}_18_1_Heo.root',\n    '${inputFilePrefixM50}_190_1_UDx.root',\n    '${inputFilePrefixM50}_191_1_S1A.root',\n    '${inputFilePrefixM50}_192_1_Bon.root',\n    '${inputFilePrefixM50}_193_1_5pX.root',\n    '${inputFilePrefixM50}_194_1_l05.root',\n    '${inputFilePrefixM50}_195_1_wyK.root',\n    '${inputFilePrefixM50}_196_1_msW.root',\n    '${inputFilePrefixM50}_197_1_5MY.root',\n    '${inputFilePrefixM50}_198_1_3qe.root',\n    '${inputFilePrefixM50}_199_1_99y.root',\n    '${inputFilePrefixM50}_19_1_GO0.root',\n    '${inputFilePrefixM50}_1_1_iX0.root',\n    '${inputFilePrefixM50}_200_1_V2h.root',\n    '${inputFilePrefixM50}_201_1_fy9.root',\n    '${inputFilePrefixM50}_202_1_NBk.root',\n    '${inputFilePrefixM50}_203_1_Rbb.root',\n    '${inputFilePrefixM50}_204_1_2CH.root',\n    '${inputFilePrefixM50}_205_1_Khn.root',\n    '${inputFilePrefixM50}_206_1_lRF.root',\n    '${inputFilePrefixM50}_207_1_CoI.root',\n    '${inputFilePrefixM50}_208_1_TP1.root',\n    '${inputFilePrefixM50}_209_1_gG6.root',\n    '${inputFilePrefixM50}_20_1_ajD.root',\n    '${inputFilePrefixM50}_210_1_klQ.root',\n    '${inputFilePrefixM50}_211_1_ibD.root',\n    '${inputFilePrefixM50}_212_1_TqO.root',\n    '${inputFilePrefixM50}_213_1_Ci5.root',\n    '${inputFilePrefixM50}_214_1_dgz.root',\n    '${inputFilePrefixM50}_215_1_EfM.root',\n    '${inputFilePrefixM50}_216_1_7ms.root',\n    '${inputFilePrefixM50}_217_1_MfK.root',\n    '${inputFilePrefixM50}_218_1_HMW.root',\n    '${inputFilePrefixM50}_219_1_dvC.root',\n    '${inputFilePrefixM50}_21_1_MNE.root',\n    '${inputFilePrefixM50}_220_1_1vu.root',\n    '${inputFilePrefixM50}_221_1_EmE.root',\n    '${inputFilePrefixM50}_222_1_7kU.root',\n    '${inputFilePrefixM50}_223_1_wY9.root',\n    '${inputFilePrefixM50}_224_1_vky.root',\n    '${inputFilePrefixM50}_225_1_kP3.root',\n    '${inputFilePrefixM50}_226_1_Qvs.root',\n    '${inputFilePrefixM50}_227_1_Xb2.root',\n    '${inputFilePrefixM50}_228_1_bbA.root',\n    '${inputFilePrefixM50}_229_2_UTg.root',\n    '${inputFilePrefixM50}_22_1_skl.root',\n    '${inputFilePrefixM50}_230_1_iDJ.root',\n    '${inputFilePrefixM50}_231_1_FG6.root',\n    '${inputFilePrefixM50}_232_1_FyJ.root',\n    '${inputFilePrefixM50}_233_1_VL1.root',\n    '${inputFilePrefixM50}_234_1_8oV.root',\n    '${inputFilePrefixM50}_235_1_2wn.root',\n    '${inputFilePrefixM50}_236_1_9oz.root',\n    '${inputFilePrefixM50}_237_1_HI4.root',\n    '${inputFilePrefixM50}_238_1_kS3.root',\n    '${inputFilePrefixM50}_239_1_9qf.root',\n    '${inputFilePrefixM50}_23_1_A3H.root',\n    '${inputFilePrefixM50}_240_1_lrR.root',\n    '${inputFilePrefixM50}_241_1_KhP.root',\n    '${inputFilePrefixM50}_242_1_3yN.root',\n    '${inputFilePrefixM50}_243_1_Yux.root',\n    '${inputFilePrefixM50}_244_1_1Jh.root',\n    '${inputFilePrefixM50}_245_1_8p7.root',\n    '${inputFilePrefixM50}_246_1_IXV.root',\n    '${inputFilePrefixM50}_247_1_t8o.root',\n    '${inputFilePrefixM50}_248_1_J17.root',\n    '${inputFilePrefixM50}_249_1_qDS.root',\n    '${inputFilePrefixM50}_24_1_BUB.root',\n    '${inputFilePrefixM50}_250_1_Pgv.root',\n    '${inputFilePrefixM50}_251_1_sd4.root',\n    '${inputFilePrefixM50}_252_1_9zW.root',\n    '${inputFilePrefixM50}_253_1_6MX.root',\n    '${inputFilePrefixM50}_254_1_Mt0.root',\n    '${inputFilePrefixM50}_255_1_7tQ.root',\n    '${inputFilePrefixM50}_256_1_vyG.root',\n    '${inputFilePrefixM50}_257_1_iev.root',\n    '${inputFilePrefixM50}_258_1_CHV.root',\n    '${inputFilePrefixM50}_259_1_pZO.root',\n    '${inputFilePrefixM50}_25_1_NlC.root',\n    '${inputFilePrefixM50}_260_1_p9n.root',\n    '${inputFilePrefixM50}_261_1_msy.root',\n    '${inputFilePrefixM50}_262_1_Tj9.root',\n    '${inputFilePrefixM50}_263_1_CSq.root',\n    '${inputFilePrefixM50}_264_1_o0R.root',\n    '${inputFilePrefixM50}_265_1_wwU.root',\n    '${inputFilePrefixM50}_266_1_RuX.root',\n    '${inputFilePrefixM50}_267_1_kjm.root',\n    '${inputFilePrefixM50}_268_1_WXP.root',\n    '${inputFilePrefixM50}_269_1_Drg.root',\n    '${inputFilePrefixM50}_26_1_67B.root',\n    '${inputFilePrefixM50}_270_1_k2y.root',\n    '${inputFilePrefixM50}_271_1_Aa8.root',\n    '${inputFilePrefixM50}_272_1_dhA.root',\n    '${inputFilePrefixM50}_273_1_NmU.root',\n    '${inputFilePrefixM50}_274_1_yw3.root',\n    '${inputFilePrefixM50}_275_1_I6p.root',\n    '${inputFilePrefixM50}_276_1_VTR.root',\n    '${inputFilePrefixM50}_277_2_6fK.root',\n    '${inputFilePrefixM50}_278_1_LvL.root',\n    '${inputFilePrefixM50}_279_1_kYq.root',\n    '${inputFilePrefixM50}_27_1_uze.root',\n    '${inputFilePrefixM50}_280_1_kwx.root',\n    '${inputFilePrefixM50}_281_1_xHL.root',\n    '${inputFilePrefixM50}_282_1_Iz7.root',\n    '${inputFilePrefixM50}_283_1_ABb.root',\n    '${inputFilePrefixM50}_284_1_Abg.root',\n    '${inputFilePrefixM50}_285_1_aWQ.root',\n    '${inputFilePrefixM50}_286_1_Ytn.root',\n    '${inputFilePrefixM50}_287_1_rUU.root',\n    '${inputFilePrefixM50}_288_1_gXs.root',\n    '${inputFilePrefixM50}_289_1_DbQ.root',\n    '${inputFilePrefixM50}_28_1_JNo.root',\n    '${inputFilePrefixM50}_290_1_pqO.root',\n    '${inputFilePrefixM50}_291_1_eNh.root',\n    '${inputFilePrefixM50}_292_1_jOZ.root',\n    '${inputFilePrefixM50}_293_1_cgH.root',\n    '${inputFilePrefixM50}_294_1_cdr.root',\n    '${inputFilePrefixM50}_295_1_MHG.root',\n    ])" "readFiles.extend([\n    '${inputFilePrefixM50}_296_1_e9L.root',\n    '${inputFilePrefixM50}_297_1_CEJ.root',\n    '${inputFilePrefixM50}_298_1_bQC.root',\n    '${inputFilePrefixM50}_299_1_jMY.root',\n    '${inputFilePrefixM50}_29_1_oxj.root',\n    '${inputFilePrefixM50}_2_1_0uT.root',\n    '${inputFilePrefixM50}_300_1_HAO.root',\n    '${inputFilePrefixM50}_301_1_cOn.root',\n    '${inputFilePrefixM50}_302_1_lj5.root',\n    '${inputFilePrefixM50}_303_1_TvW.root',\n    '${inputFilePrefixM50}_304_1_4Ra.root',\n    '${inputFilePrefixM50}_305_1_6mS.root',\n    '${inputFilePrefixM50}_306_1_deL.root',\n    '${inputFilePrefixM50}_307_1_THW.root',\n    '${inputFilePrefixM50}_308_1_8HZ.root',\n    '${inputFilePrefixM50}_309_1_XTF.root',\n    '${inputFilePrefixM50}_30_1_wly.root',\n    '${inputFilePrefixM50}_31_1_doP.root',\n    '${inputFilePrefixM50}_32_1_H7w.root',\n    '${inputFilePrefixM50}_33_1_Bax.root',\n    '${inputFilePrefixM50}_34_1_Rbe.root',\n    '${inputFilePrefixM50}_35_1_Jc4.root',\n    '${inputFilePrefixM50}_36_1_nNT.root',\n    '${inputFilePrefixM50}_37_1_J8R.root',\n    '${inputFilePrefixM50}_38_1_Mqf.root',\n    '${inputFilePrefixM50}_39_1_pc0.root',\n    '${inputFilePrefixM50}_3_1_GFw.root',\n    '${inputFilePrefixM50}_40_1_gQ0.root',\n    '${inputFilePrefixM50}_41_1_llv.root',\n    '${inputFilePrefixM50}_42_1_DGn.root',\n    '${inputFilePrefixM50}_43_1_Nby.root',\n    '${inputFilePrefixM50}_44_1_7E2.root',\n    '${inputFilePrefixM50}_45_1_34W.root',\n    '${inputFilePrefixM50}_46_1_hrg.root',\n    '${inputFilePrefixM50}_47_1_Y2k.root',\n    '${inputFilePrefixM50}_48_1_bIV.root',\n    '${inputFilePrefixM50}_49_1_plU.root',\n    '${inputFilePrefixM50}_4_1_GYc.root',\n    '${inputFilePrefixM50}_50_1_jdg.root',\n    '${inputFilePrefixM50}_51_1_RH4.root',\n    '${inputFilePrefixM50}_52_1_f7c.root',\n    '${inputFilePrefixM50}_53_1_neH.root',\n    '${inputFilePrefixM50}_54_1_gte.root',\n    '${inputFilePrefixM50}_55_1_u97.root',\n    '${inputFilePrefixM50}_56_1_zve.root',\n    '${inputFilePrefixM50}_57_1_bTz.root',\n    '${inputFilePrefixM50}_58_1_MBq.root',\n    '${inputFilePrefixM50}_59_3_CbD.root',\n    '${inputFilePrefixM50}_5_1_QO5.root',\n    '${inputFilePrefixM50}_60_1_PA2.root',\n    '${inputFilePrefixM50}_61_1_RGL.root',\n    '${inputFilePrefixM50}_62_1_Coj.root',\n    '${inputFilePrefixM50}_63_1_67s.root',\n    '${inputFilePrefixM50}_64_1_ksM.root',\n    '${inputFilePrefixM50}_65_1_nEi.root',\n    '${inputFilePrefixM50}_66_1_SnN.root',\n    '${inputFilePrefixM50}_67_1_HkD.root',\n    '${inputFilePrefixM50}_68_1_Fso.root',\n    '${inputFilePrefixM50}_69_1_QJr.root',\n    '${inputFilePrefixM50}_6_1_D9t.root',\n    '${inputFilePrefixM50}_70_1_LbP.root',\n    '${inputFilePrefixM50}_71_1_YLo.root',\n    '${inputFilePrefixM50}_72_1_EAY.root',\n    '${inputFilePrefixM50}_73_1_Gwy.root',\n    '${inputFilePrefixM50}_74_1_PLq.root',\n    '${inputFilePrefixM50}_75_1_bJ8.root',\n    '${inputFilePrefixM50}_76_1_GhV.root',\n    '${inputFilePrefixM50}_77_1_W8h.root',\n    '${inputFilePrefixM50}_78_1_ZKE.root',\n    '${inputFilePrefixM50}_79_1_ZMW.root',\n    '${inputFilePrefixM50}_7_1_rDM.root',\n    '${inputFilePrefixM50}_80_1_CiJ.root',\n    '${inputFilePrefixM50}_81_1_HXw.root',\n    '${inputFilePrefixM50}_82_1_eJO.root',\n    '${inputFilePrefixM50}_83_1_G5h.root',\n    '${inputFilePrefixM50}_84_1_l1O.root',\n    '${inputFilePrefixM50}_85_1_Y8b.root',\n    '${inputFilePrefixM50}_86_1_fEs.root',\n    '${inputFilePrefixM50}_87_1_0YW.root',\n    '${inputFilePrefixM50}_88_1_4i2.root',\n    '${inputFilePrefixM50}_89_1_W0B.root',\n    '${inputFilePrefixM50}_8_1_XrF.root',\n    '${inputFilePrefixM50}_90_1_dgj.root',\n    '${inputFilePrefixM50}_91_1_Knd.root',\n    '${inputFilePrefixM50}_92_1_hQz.root',\n    '${inputFilePrefixM50}_93_1_gMI.root',\n    '${inputFilePrefixM50}_94_1_Nu8.root',\n    '${inputFilePrefixM50}_95_1_P9c.root',\n    '${inputFilePrefixM50}_96_1_2QX.root',\n    '${inputFilePrefixM50}_97_1_UMG.root',\n    '${inputFilePrefixM50}_98_1_Dp9.root',\n    '${inputFilePrefixM50}_99_1_0jW.root',\n    '${inputFilePrefixM50}_9_1_sde.root',\n    ])" )

#CleanJets output file
cleanJetsOutFiles=( "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_NonIsoWDYJetsToLL_M-10To50.root" "${cleanJetsOutputFilePrefix}NMSSMSignal_MuProperties_NonIsoWDYJetsToLL_M-50.root" )

#TauAnalyzer output files
highMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_NonIsoWDYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_highMT_NonIsoWDYJetsToLL_M-50_${version}.root" )
lowMTIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_NonIsoWDYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadIsoAnalysis_lowMT_NonIsoWDYJetsToLL_M-50_${version}.root" )
highMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_NonIsoWDYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_highMT_NonIsoWDYJetsToLL_M-50_${version}.root" )
lowMTNonIsoTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_NonIsoWDYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadNonIsoAnalysis_lowMT_NonIsoWDYJetsToLL_M-50_${version}.root" )
highMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_NonIsoWDYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_highMT_NonIsoWDYJetsToLL_M-50_${version}.root" )
lowMTAllTauAnalyzerOutputFiles=( "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_NonIsoWDYJetsToLL_M-10To50_${version}.root" "${tauAnalyzerOutputFilePrefix}muHadAnalysis_lowMT_NonIsoWDYJetsToLL_M-50_${version}.root" )

#EDM output files
EDMOutputFiles=( "${EDMOutputFilePrefix}NonIsoWDYJetsToLL_M-10To50${infoTag}_${version}.root" "${EDMOutputFilePrefix}NonIsoWDYJetsToLL_M-50${infoTag}_${version}.root" )

#samples
samples=( "NonIsoWDYJetsToLL_M-10To50" "NonIsoWDYJetsToLL_M-50" )

####GENERATION LOOP####

#change to working directory
mkdir -p $dir
cd $dir

#loop over number of samples
for i in `seq $iBeg $iEnd`
  do

  #generate cfg file
  sed -e "s%FILES%${inputFileBlocks[${i}]}%" -e "s%CLEANJETSOUTFILE%${cleanJetsOutFiles[${i}]}%" -e "s%HIGHMTNONISOTAUANALYZEROUTFILE%${highMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTALLTAUANALYZEROUTFILE%${highMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%HIGHMTISOTAUANALYZEROUTFILE%${highMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTNONISOTAUANALYZEROUTFILE%${lowMTNonIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTALLTAUANALYZEROUTFILE%${lowMTAllTauAnalyzerOutputFiles[${i}]}%" -e "s%LOWMTISOTAUANALYZEROUTFILE%${lowMTIsoTauAnalyzerOutputFiles[${i}]}%" -e "s%EDMOUTFILE%${EDMOutputFiles[${i}]}%" -e "s%HIGGSREW%False%"  e "s%REWEIGHT%False%" -e "s%PUSCENARIO%S10%" -e "s%SAMPLE%%" ../${templateCfg} > tauanalyzer_${samples[${i}]}_cfg.py

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
cat <<EOF > runNonIsoWDYJetsToLLTauAnalyzerCfgs.sh
#!/bin/bash

for file in \`ls -alh *NonIsoWDYJetsToLL*.py | awk '{ print \$9 }'\`
  do
  outFile=\`echo \$file | sed -e "s%\.py%.txt%"\`
  cmsRun \$file > \$outFile
done

exit 0
EOF
chmod a+x runNonIsoWDYJetsToLLTauAnalyzerCfgs.sh

#generate script that submits all jobs to LSF
cat <<EOF > submitNonIsoWDYJetsToLLTauAnalyzerJobs.sh
#!/bin/bash

for file in \`ls -alh tauanalyzer*NonIsoWDYJetsToLL*.sh | awk '{ print \$9 }'\`
  do
  jobName=\`echo \$file | sed -e "s%\(.*\)\.sh%\1%"\`
  bsub -q 1nd -J \$jobName < \$file
done

exit 0
EOF
chmod a+x submitNonIsoWDYJetsToLLTauAnalyzerJobs.sh

#generate script that copies all files locally from EOS
cat <<EOF > copyNonIsoWDYJetsToLLFromEOS.sh
#!/bin/bash

eval \`scramv1 runtime -sh\`
for sample in "10To50" "50"
  do
  for cut in "Iso" "NonIso" ""
    do
    for MTBin in high low
      do
      cmsStage -f /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_NonIsoWDYJetsToLL_M-\${sample}_${version}.root /data1/`whoami`/nonIsoWDYJetsToLL/analysis/
      cmsRm /store/user/`whoami`/muHad\${cut}Analysis_\${MTBin}MT_NonIsoWDYJetsToLL_M-\${sample}_${version}.root
    done
  done
done

exit 0
EOF
chmod a+x copyNonIsoWDYJetsToLLFromEOS.sh

exit 0
