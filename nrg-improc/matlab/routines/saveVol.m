function []=saveVol(V,fname)
V32=single(V);
sz=size(V);
avw=avw_hdr_make();
avw.hdr.dime.dim(1:8) = int16([4 sz(1) sz(2) sz(3) 1 0 0 0]);
avw.hdr.dime.datatype = int16(16);
avw.hdr.dime.bitpix = int16(32);
avw.img=V32;
avw_img_write(avw,fname);