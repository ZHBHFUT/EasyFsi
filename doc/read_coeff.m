clear
clc

fid=fopen('coeff.txt','rb');

nbs=fread(fid,1,"int");
nbt=fread(fid,1,"int");

nnode_s=zeros(nbs,1);
nface_s=zeros(nbs,1);
nnode_t=zeros(nbt,1);
nface_t=zeros(nbt,1);
for k=1:nbs
    nnode_s(k)=fread(fid,1,"int");
    nface_s(k)=fread(fid,1,"int");
end
for k=1:nbt
    nnode_t(k)=fread(fid,1,"int");
    nface_t(k)=fread(fid,1,"int");
end

nn = sum(nnode_t);
nf = sum(nface_t);
max_donor=20;

node_src_bd=zeros(nn,1);
node_ndonor=zeros(nn,1);
node_donors=zeros(nn,max_donor);
node_weight=zeros(nn,max_donor);
node_drng=zeros(nn,2);
node_dist=zeros(nn,1);

for k=1:nn
    node_src_bd(k)=fread(fid,1,'int');
    node_ndonor(k)=fread(fid,1,'int');
    node_donors(k,:)=fread(fid,max_donor,'int');
    node_weight(k,:)=fread(fid,max_donor,'double');
    node_drng(k,1:2)=fread(fid,2,'uint64');
    node_dist(k)=fread(fid,1,'double');

    node_donors(k,1:node_ndonor(k))=node_donors(k,1:node_ndonor(k))+1;
end

face_src_bd=zeros(nf,1);
face_ndonor=zeros(nf,1);
face_donors=zeros(nf,max_donor);
face_weight=zeros(nf,max_donor);
face_drng=zeros(nf,2);
face_dist=zeros(nf,1);
for k=1:nf
    face_src_bd(k)=fread(fid,1,'int');
    face_ndonor(k)=fread(fid,1,'int');
    face_donors(k,:)=fread(fid,max_donor,'int');
    face_weight(k,:)=fread(fid,max_donor,'double');
    face_drng(k,1:2)=fread(fid,2,'uint64');
    face_dist(k)=fread(fid,1,'double');

    face_donors(k,1:face_ndonor(k))=face_donors(k,1:face_ndonor(k))+1;
end

fclose(fid);