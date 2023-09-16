% Comments:	Mikhail Milchenko, mmilchenko@wustl.edu, Washington University School of Medicine. 
% Algorithm:	Facemasking 1

% Mask: binary mask of VOI
% thickness: level depth
% step: partition step
% thresh: mask threshold
% ac: depth anisotropy coefficient (ratio of upper thickness to lower thickness)
% nm: strel size in pix/mm for morphology operations (~2pix for 1x1x1 mm)

function[faces_o, vertices_o, faces_b, vertices_b, faces_t,vertices_t,...
    lower_surf,upper_surf,orig_surf]=RectangularMesh(root,Mask,thickness,step,thresh,ac,nm,opt)

LEFT=1;TOP=2;RIGHT=4;BOTTOM=8;

sz=size(Mask);

%morphologically close.
se=strel('disk',nm);
Mask2=imdilate(Mask,se);

%step=thickness*2;
dx=max(1,floor((sz(1)-1)/step(1)));
dy=max(1,floor((sz(2)-1)/step(2)));

%generate faces
for x=1:dx
    codex=0;
    lx=(x-1)*step(1)+1;
    if(x==1) codex=LEFT;
    elseif(x==dx) codex=RIGHT;
    else codex=0; end;
    for y=1:dy
        codey=0;
        if(y==1) codey=TOP;
        elseif(y==dy) codey=BOTTOM;
        else codey=0;
        end;
        ty=(y-1)*step(2)+1;
        faces(x,y)=qface([lx,ty],[step(1), step(2)],sz,Mask,thresh,[x,y],bitor(codex,codey));
    end;
end;

for x=1:dx+1
    if(x==1) lt=[]; lb=[]; end;
    if(x==dx+1) rt=[]; rb=[]; end;        
    for y=1:dy+1
        if(y==1) lt=[]; rt=[]; end;
        if(y==dy+1) lb=[]; rb=[]; end;
        
        if(x>1 && y>1) lt=faces(x-1,y-1);end;
        if(x<dx+1 && y>1) rt=faces(x,y-1);end;
        if(x>1 && y<dy+1) lb=faces(x-1,y);end;
        if(x<dx+1 && y<dy+1) rb=faces(x,y);end;
        vertices(x,y,:)=vertex(lt,rt,lb,rb,thickness,ac);
    end;    
end;

%improve coverage of the upper surface.
[vertices Mask1]=inflate_surf(vertices,dx,dy,Mask2,opt);
if(opt>1)
    saveVol(double(Mask1),[root '_intersect_layer']);
    Mask1=raster_layer(vertices,dx,dy,Mask2);
    saveVol(double(Mask1),[root '_raster_layer']);
end;

%C=calc_curvature(vertices,dx,dy);
lower_surf=[];
upper_surf=[];
orig_surf=[];
for x=1:dx+1
    lx=(x-1)*step(1)+1;
    for y=1:dy+1
        ty=(y-1)*step(2)+1;
        if(vertices(x,y,:).empty()~=1)
          lower_surf(:,x,y)=vertices(x,y).xb;
          upper_surf(:,x,y)=vertices(x,y).xt;
          orig_surf(:,x,y)=vertices(x,y).x;
        else
            lower_surf(:,x,y)=[lx,ty,1];
            upper_surf(:,x,y)=[lx,ty,1];
            orig_surf(:,x,y)=[lx,ty,1];
        end;
    end;
end;

[faces_b,vertices_b]=generate_mesh(lower_surf,orig_surf);
[faces_t,vertices_t]=generate_mesh(upper_surf,orig_surf);
[faces_o,vertices_o]=generate_mesh(orig_surf,orig_surf);

return;
%end of RectangularMesh function


function [C] = calc_curvature(vertices,dx,dy)
    for x=1:dx
        for y=1:dy
            code=0;
            n=0;
            k=0;
            v0=vertices(x,y);
            if(x==1) code=bitor(code,1);end;
            if(x==dx) code=bitor(code,2);end;
            if(y==1) code=bitor(code,4);end;
            if(y==dy) code=bitor(code,8);end;
            if(bitand(code,1)==0 && bitand(code,2)==0)
                k=max(k,acos(ang(v0,vertices(x-1,y),vertices(x+1,y))));
            end;
            if(bitand(code,4)==0 && bitand(code,8)==0)
                k=max(k,acos(ang(v0,vertices(x,y-1),vertices(x,y+1))));
            end;
            if(bitand(code,4)==0 && bitand(code,2)==0)
                k=max(k,acos(ang(v0,vertices(x,y-1),vertices(x+1,y))));
            end;
            if(bitand(code,8)==0 && bitand(code,2)==0)
                k=max(k,acos(ang(v0,vertices(x+1,y),vertices(x,y+1))));
            end;
            if(bitand(code,8)==0 && bitand(code,1)==0)
                k=max(k,acos(ang(v0,vertices(x-1,y),vertices(x,y+1))));
            end;
            if(bitand(code,1)==0 && bitand(code,4)==0)
                k=max(k,acos(ang(v0,vertices(x-1,y),vertices(x,y-1))));
            end;                        
            C(x,y)=k;
        end;
    end;
return;

function [c]=curv(v0,v1)
n=v1.xb-v0.xb;
nm=norm(n);
if(nm>0) c=1/nm;
else c=0; end;
return;

function [a]=ang(v0,v1,v2)
a1=v1.xb-v0.xb; a2=v2.xb-v0.xb;
if(norm(a1)==0 || norm(a2)==0) a=1;
else a=abs(dot(a1,a2)/(norm(a1)*norm(a2)));
end;
return;

function [Mask1]=raster_layer(vertices,dx,dy,Mask)
Mask1=zeros(size(Mask));
Maskz=zeros(size(Mask));
for x=1:dx
    for y=1:dy
        if(vertices(x,y,:).empty()~=1)
          vert0=vertices(x,y);
          vertr=vertices(x+1,y);
          vertb=vertices(x,y+1);
          vertbr=vertices(x+1,y+1);
          [r Mask1]=triatest(vert0.xt,vertr.xt,vertbr.xt,Maskz,Mask1,4);
          [r Mask1]=triatest(vert0.x,vertr.x,vertbr.x,Maskz,Mask1,3);
          [r Mask1]=triatest(vert0.xb,vertr.xb,vertbr.xb,Maskz,Mask1,2);
          
          [r Mask1]=triatest(vert0.xt,vertb.xt,vertbr.xt,Maskz,Mask1,4);
          [r Mask1]=triatest(vert0.x,vertb.x,vertbr.x,Maskz,Mask1,3);
          [r Mask1]=triatest(vert0.xb,vertb.xb,vertbr.xb,Maskz,Mask1,2);
        end;
    end;
end;
return;

function[new_vert Mask1]=inflate_surf(vertices,dx,dy,Mask,opt)
if(opt==0)
    Mask1=[];
else
    Mask1=zeros(size(Mask));
end;

for niter=1:10
    nmod=0;
    for x=1:dx
        for y=1:dy
            if(vertices(x,y,:).empty()~=1)
                vert0=vertices(x,y);
                vertr=vertices(x+1,y);
                vertb=vertices(x,y+1);
                vertbr=vertices(x+1,y+1);
                [r Mask1]=triatest(vert0.xt,vertr.xt,vertbr.xt,Mask,Mask1,opt);
                if(r>0)
                    vertices(x,y).mod=1; %1
                    vertices(x+1,y).mod=1; %2
                    vertices(x+1,y+1).mod=1; %3                   
                    nmod=nmod+1;
                end;
                if(opt>1)
                    [r Mask1]=triatest(vert0.xb,vertr.xb,vertbr.xb,Mask,Mask1,opt);
                end;
                [r Mask1]=triatest(vert0.xt,vertb.xt,vertbr.xt,Mask,Mask1,opt);
                if(r>0)
                    vertices(x,y).mod=1; %1
                    vertices(x,y+1).mod=1; %4
                    vertices(x+1,y+1).mod=1; %3                   
                    nmod=nmod+1;
                end;
                if(opt>1)
                    [r Mask1]=triatest(vert0.xb,vertb.xb,vertbr.xb,Mask,Mask1,opt);
                end;
            end;
        end;
    end;
    if(nmod==0) break; end;
    %now update all vertices
    for x=1:dx+1
        for y=1:dy+1
            if (vertices(x,y).mod>0)
                temp=vertices(x,y);
                vertices(x,y)=stepup(temp);
                vertices(x,y).mod=0;
            end;
        end;
   end;
end;
new_vert=vertices;
return;

function[res,Mask1]=triatest(x0,x1,x2,Mask,Mask1,opt)
v1=x1-x0;
v2=x2-x0;
ndiv1=ceil(abs(max(v1(:))));
ndiv2=ceil(abs(max(v2(:))));
ndiv=max(ndiv1,ndiv2);
dv1=v1/ndiv;
dv2=v2/ndiv;
res=0;
res=max(res,testinter(x0,x1,Mask,Mask1,opt));
if(opt<2 && res>0) return; end;
res=max(res,testinter(x0,x2,Mask,Mask1,opt));
if(opt<2 && res>0) return; end;
res=max(res,testinter(x1,x2,Mask,Mask1,opt));
if(opt<2 && res>0) return; end;

for l=0:ndiv
    x_1=x0+dv1*l;
    x_2=x0+dv2*l;
    [res1 Mask1]=testinter(x_1,x_2,Mask,Mask1,opt);
    res=max(res,res1);
    if(opt<2 && res>0) return; end;
end;
    
% test for intersections of an edge with mask.
function[res Mask1]=testinter(x0,x1,Mask,Mask1,opt)
sz=int16(size(Mask));
v=x1-x0;
ndiv=ceil(abs(max(v(:))));
dv=v/ndiv;
res=0;
%Mask1=[];
for l=0:ndiv
    x=int16(round(x0+dv*l));
    if(sum(x<=0)>0) continue; end;
    if(sum(x>sz)>0) continue; end;
    if(opt>1)
        Mask1(x(1),x(2),x(3))=opt;
    end;
    if (Mask(x(1),x(2),x(3)) ~= 0)        
        res=1;
        if(opt>1)
            Mask1(x(1),x(2),x(3))=opt;    
        end;
        if(opt<2) return; end;
    end
end
return;

function [vret]=stepup(vert)
vret=vert;
if (vert.nmod>=10) return; end;
vret.mod=1;
vret.nmod=vret.nmod+1;
dv=norm(vret.x-vret.xb);
vret.xt=vret.xt+vret.n*dv*.1;
return;

function[new_pt]=modify_surf_pt(upper_surf, reg_surf, normals,x,y)
max_r=2.0;
pt0=reg_surf(:,x,y);
pleft=upper_surf(:,x-1,y);
ptop=upper_surf(:,x,y-1);
pright=upper_surf(:,x+1,y);
pbottom=upper_surf(:,x,y+1);
pt=upper_surf(:,x,y);
plane1=plane_eq(pright,pleft,ptop);
plane2=plane_eq(pleft, pright,pbottom);
res1=0;res2=0;
if(dot([pt(1),pt(2),pt(3),1],plane1)<-1e-6) res1=1;
else res1=0; end;
if(dot([pt(1),pt(2),pt(3),1],plane2)<-1e-6) res2=1;
else res2=0; end;


if(max(plane1)>1e-3 && res1~=0)    
    [pt1,res1]=intersect_line_with_plane(pt0,normals(:,x,y),plane1);
end;
if(max(plane2)>1e-3 && res2~=0)
    [pt2,res2]=intersect_line_with_plane(pt0,normals(:,x,y),plane2);
end;
if(res1==0 && res2==0) new_pt=pt;return;end;

n=norm(pt-pt0);
if(n<1e-6) r0=0; n=1; 
else r0=1;end;

if(res1~=0) r1=norm(pt1'-pt0)/n;
else r1=0; end;
if(res2~=0) r2=norm(pt2'-pt0)/n;
else r2=0; end;

if(r1>max_r) r1=0;end;
if(r2>max_r) r2=0;end;

if(r1>r2 && r1>0) new_pt=pt1';
elseif(r2>=r1 && r2>0) new_pt=pt2';
else new_pt=pt; end;
return;

function[new_pt,res]=intersect_line_with_plane(pt0,normal,c)
cx=c(1);cy=c(2);cz=c(3);c0=c(4);
nx=normal(1);ny=normal(2);nz=normal(3);
den=cx*nx+cy*ny+cz*nz;
if(abs(den)<1e-6) res=0;new_pt=[0,0,0];return;end;
x0=pt0(1);y0=pt0(2);z0=pt0(3);
x=-(-ny*x0*cy-x0*cz*nz+nx*cz*z0+nx*c0+nx*cy*y0)/den;
y=-(cx*x0*ny-cx*nx*y0-cz*nz*y0+cz*z0*ny+c0*ny)/den;
z=(-nz*cx*x0-nz*c0-nz*cy*y0+z0*cx*nx+ny*z0*cy)/den;
res=1;
new_pt=[x,y,z];
return;

function[faces,vertcs]=generate_mesh(surf,orig_surf)
sz=size(surf);
dx=sz(2)-1; dy=sz(3)-1;

face_ind=1;
vert_ind=1;
faces=[]; 
vertcs=[];
for y=1:dy
    for x=1:dx
        
        lt=surf(:,x,y);
        ltb=orig_surf(:,x,y);
        
        rt=surf(:,x+1,y);
        rtb=orig_surf(:,x+1,y);
        
        rb=surf(:,x+1,y+1);
        rbb=orig_surf(:,x+1,y+1);
        
        lb=surf(:,x,y+1);
        lbb=orig_surf(:,x,y+1);
        
        [vrt,fcs]=qpatch(lt,rt,rb,lb,ltb,rtb,rbb,lbb);
%        [vrtb,fcsb]=qpatch(ltb,rtb,rbb,lbb);
        
        if(isempty(vrt)) continue;end;
        
        nvert=size(vrt,1);
        nface=size(fcs,1);
        vertcs(:,vert_ind:vert_ind+nvert-1)=vrt';
        faces(face_ind,:)=(fcs(1,:)+vert_ind-1);
        if(nface==2)
            faces(face_ind+1,:)=fcs(2,:)+vert_ind-1;
        end;
        vert_ind=vert_ind+nvert;
        face_ind=face_ind+nface;
        
%         if(surf(3,x,y)>1 || surf(3,x,y+1)>1 || surf(3,x+1,y)>1 || surf(3,x+1,y+1)>1)
%             vertices(:,vert_ind)=surf(:,x,y);
%             vertices(:,vert_ind+1)=surf(:,x,y+1);
%             vertices(:,vert_ind+2)=surf(:,x+1,y+1);
%             vertices(:,vert_ind+3)=surf(:,x+1,y);
%             faces(:,face_ind)=[vert_ind,vert_ind+1,vert_ind+2,vert_ind+3];
%             face_ind=face_ind+1;
%             vert_ind=vert_ind+4;
%         end;

    end;
end;
faces=faces';
return;

function[vertcs,faces]=qpatch(lt,rt,rb,lb,ltb,rtb,rbb,lbb)
nvert=0;
if(ltb(3)~=1)
    nvert=nvert+1;        
    vertcs(nvert,:)=lt;
end;
if(rtb(3)~=1)
    nvert=nvert+1;
    vertcs(nvert,:)=rt;
end;
if(rbb(3)~=1)
    nvert=nvert+1;    
    vertcs(nvert,:)=rb;
end;
if(lbb(3)~=1)
    nvert=nvert+1;
    vertcs(nvert,:)=lb;
end;
if(nvert==3)
    faces=[1,2,3];
    return;
elseif(nvert==4)
    faces=[[1,2,3];[3,4,1]];
else
    faces=[];
    vertcs=[];
    return;
end;

function[dest_surf]=filter_surf(src_surf,model_surf,step,thick)
dest_surf=src_surf;
sz=size(src_surf);
dx=sz(2);dy=sz(3);
for y=1:dy
    for x=1:dx
        dest_surf(:,x,y)=allowable_edge(src_surf(:,x,y),model_surf(:,x,y),...
            step,thick,x,y);
    end;
end;
return;

function[pt]=allowable_edge(pt0,pt1,step,thick,x,y)
if(norm(pt0-pt1)<max(step,thick*2+1)*5) pt=pt0;
else
    pt=[(x-1)*step+1,(y-1)*step+1,1];
end;
return;

%averaged unit normal to surface in point (x,y).
function[n]=calc_avg_normal(normals, x,dx,y,dy)
n1=normals(:,1,x,y);

if(y>1)
    n2=normals(:,2,x,y-1);
else
    n2=[0 0 1]';
end;

if(y>1 && x>1)
    n3=normals(:,3,x-1,y-1);
else
    n3=[0 0 1]';
end;

if(x>1)
    n4=normals(:,4,x-1,y);
else
    n4=[0 0 1]';
end;
n=(n1+n2+n3+n4);
n=n/norm(n);

%return 0 if vector is not within bounds; 1 otherwise
function[res]=ib(v,sz, nD)
res=1;
for i=1:nD
    if(v(i)<1) res=0; return;
    elseif(v(i)>sz(i)) res=0; return; end;
end;
%end of ib (iN bOUNDS) function.
    
%Unit normal to plane defined by three points.
function[x]=unit_normal(x0,x1,x2)
x=cross(x0-x1,x2-x1);
x=x/norm(x);
function [f] = face(center,p1,p2)
f.center=center;
f.p1=p1;
f.p2=p2;
n=cross(p1-center,p2-center);
f.n=n/norm(n);


function [xt,xb] = calc_edge(faces,h)
nfaces=size(faces,1);
n=faces(1).n;
for i=2:nfaces
    n=n+faces(i).n;
end;

n=n/norm(n);
den=0;
for i=1:nfaces
    den=den+1/dot(n,faces(i).n);
end;
vedge=(1.0/nfaces)*h*den*n;   
xt=faces(1).center+vedge;
xb=faces(1).center-vedge;

%using pre-calculated points of face, return adjusted points in cross vicinity.
function [faces] = calc_faces(x,y,dx,dy,step,tl,bl,tr,br)

lx=(x-2)*step+1;
mx=(x-1)*step+1;
rx=x*step+1;

ty=(y-2)*step+1;
my=(y-1)*step+1;
by=y*step+1;

ind=1;

if(x<dx && y<dy)
   faces(ind)=face([mx,my,tl(x,y)],[rx,ty,tr(x,y)],[mx,by,br(x,y)]);
   ind=ind+1;
end
if(x>1 && y<dy)
   faces(ind)=face([mx,my,tl(x,y)],[mx,by,br(x,y)],[lx,my,tl(x-1,y)]);
   ind=ind+1;
end
if(x>1 && y>1)
    faces(ind)=face([mx,my,tl(x,y)],[lx,my,tl(x-1,y)],[mx,ty,tl(x,y-1)]);
    ind=ind+1;
end;
if(x<dx && y>1)
    faces(ind)=face([mx,my,tl(x,y)],[mx,ty,tl(x,y-1)],[rx,my,tr(x,y)]);
    ind=ind+1;
end;


%using pre-calculated points of face, return adjusted points of the same grid element.
function [ptl,pbl,ptr,pbr] = calc_pt_sq(x,y,dx,dy,step,tl,bl,tr,br)
lx=(x-1)*step+1;
rx=lx+step;
ty=(y-1)*step+1;
by=ty+step;
if(x<=dx && y<=dy)
    vtl=tl(x,y);
    vbl=bl(x,y);
    vtr=tr(x,y);
    vbr=br(x,y);
elseif(x==dx+1 && y==dy+1)
    vtl=br(x-1,y-1);
    vbl=1; vtr=1;vbr=1;
elseif(x==dx+1)
    vtl=tr(x-1,y); vbl=br(x-1,y);
    vtr=1; vbr=1;
elseif(y==dy+1)
    vtl=bl(x,y-1); vbl=1;
    vtr=br(x,y-1); vbr=1;
else
    vtl=1;vbl=1;vtr=1;vbr=1;
end;
ptl=[lx,ty,vtl]; pbl=[lx,by,vbl];
ptr=[rx,ty,vtr]; pbr=[rx,by,vbr];

%v0 is lower end of edge through x1
%v1 is upper end of edge through x1
%x0,x1,x2 - face vertices
%v - averaged normal at x1
%w - layer thickness
%sz - volume bounds
function[v0,v1]=edge_vector(x0,x1,x2,v,w,sz)
if(x1(3)==1)
    v0=x1;
    v1=x1;
    return;
end;
n=unit_normal(x0,x1,x2);
cosf=dot(n,v);
n=v*w/cosf;
v0=x1-n';
v1=x1+n';
%v0=atbv(x1-n,sz);
%v1=atbv(x1+n,sz);

%performs rectangular ray cast from maximum z coordinate to
%minimum z coordinate, searching for a non-zero value.

function[res]=getval(Mask,x,axis,sz,thresh)
if(axis>0) dir=1; else dir=-1; axis=-axis; end;

if(dir==-1) st=sz(axis);en=1;res=1;
else st=1; en=sz(axis);res=sz(axis);
end;

for i=st:dir:en
    switch(axis)
        case 3
            val=Mask(x(1),x(2),i);
        case 2 
            val=Mask(x(1),i,x(3));
        case 1 
            val=Mask(i,x(2),x(3));
    end;            
    if(val>thresh)
      res=i;
      break;
    end;
end;
return;

function[vert] = vertex(lt,rt,lb,rb,h,ac)
nNorm=0;
% initialize vert (should be done here).
vert.x=[0,0,0];
vert.mod=0;
vert.nmod=0;
vert.empty=1;
vert.n=[0,0,0];
vert.xt=[0,0,0];
vert.xb=[0,0,0];

if(~isempty(lt))
    if(lt.base==0)
        vert.x=lt.rb;
        normals(nNorm+1,:)=lt.nrb;
        nNorm=nNorm+1;
    end;
end;
if(~isempty(rt))
    if(rt.base==0)
        vert.x=rt.lb;
        normals(nNorm+1,:)=rt.nlb;
        nNorm=nNorm+1;
    end;
end;
if(~isempty(rb))
    if(rb.base==0)
        vert.x=rb.lt;
        normals(nNorm+1,:)=rb.nlt;
        nNorm=nNorm+1;
    end;
end;
if(~isempty(lb))
    if(lb.base==0)
        vert.x=lb.rt;
        normals(nNorm+1,:)=lb.nrt;
        nNorm=nNorm+1;
    end;
end;
vert.mod=0;
vert.nmod=0;
if(nNorm==0) return;
else vert.empty=0;end;

n=normals(1,:);
for(i=2:nNorm)
    n=n+normals(i,:);
end;
n=n/norm(n);
vert.n=n;

sum=0;
for(i=1:nNorm)
    sum=sum+1/dot(n,normals(i,:));
end;
%n=n*(1/nNorm)*sum*h;
n=n*h;
vert.xt=vert.x+ac*n+(rand()-0.5)*8;
vert.xb=vert.x-n+(rand()-0.5)*8;

function[vert] = fvertex(x,y,codex,codey,Mask,sz,thresh)
vert.x=[l,t,getval(Mask,[x,y],-3,sz,thresh)];
vx=getval(Mask,[x,y],codex,sz,thresh);
vy=getval(Mask,[x,y],codey,sz,thresh);
if(vx==1 || vy==1) vert.base=1; else vert.base=0;end;
    

function[face] = qface(x,dx,sz,Mask,thresh,ind,code)

l=x(1);r=x(1)+dx(1);
t=x(2);b=x(2)+dx(2);
face.code=code;
face.ix=ind(1);
face.iy=ind(2);

face.lt=[l,t,getval(Mask,[l,t],-3,sz,thresh)];
face.rt=[r,t,getval(Mask,[r,t],-3,sz,thresh)];
face.lb=[l,b,getval(Mask,[l,b],-3,sz,thresh)];
face.rb=[r,b,getval(Mask,[r,b],-3,sz,thresh)];

if(face.lt(3)==face.rt(3)==face.lb(3)==face.rb(3)==1)
   face.base=1;
elseif(face.lt(3)==1 || face.rt(3)==1 || face.lb(3)==1 || face.rb(3)==1)
   face.base=0;
else
    face.base=0;
end;
face.nlt=unit_normal(face.rt,face.lt,face.lb);
face.nrt=unit_normal(face.rb,face.rt,face.lt);
face.nrb=unit_normal(face.lb,face.rb,face.rt);
face.nlb=unit_normal(face.lt,face.lb,face.rb);
