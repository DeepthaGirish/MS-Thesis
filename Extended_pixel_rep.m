 
clc;
 clear all;
close all;
% Create equation for number of clusters
% data='imagedata.xlsx';
% x= xlsread(data,'A:C');
% y=xlsread(data,'D:D');
%w=[5.6955 16];
% w=inv(transpose(x)*x)*(transpose(x))*y;
x1=imread('D:\University of Cincinnati\MS Thesis\EPR project\EPR data\test.jpg');
 
x= imresize(x1,.5);
x=im2double(x);
xhsv=rgb2hsv(x1);
xhsv=imresize(xhsv,1);
xhsv=im2double(xhsv);
xy=rgb2ycbcr(x1);
xy=imresize(xy,1);
xy=im2double(xy);
big=size(x);
row=big(1);
col=big(2);
y=zeros(row,col,3);
eprx(:,:,:)=0;
clusters(:,:,:)=0;
% clno=w(1)+(w(2).*std2(x));
% clno=round(clno);
% clno=clno+2;
clno=5;
 
 
clusters(:,:,:);
mclu(:,:)=0;% holds mean vector for every cluster
R=xy(:,:,1);G=xy(:,:,2);B=xy(:,:,3);% eprx is the new pixel rep as described in mam's paper
 
for i=1:1:row
    for j=1:1:col
        eprx(i,j,1)=R(i,j);
        eprx(i,j,2)=G(i,j);
        eprx(i,j,3)=B(i,j);
        eprx(i,j,4)=sqrt((((R(i,j))^2)+((G(i,j))^2)));
         eprx(i,j,5)=atan((G(i,j)/R(i,j)));
        eprx(i,j,6)=sqrt((((G(i,j))^2)+((B(i,j))^2)));
         eprx(i,j,7)=atan((B(i,j)/G(i,j)));
        eprx(i,j,8)=sqrt((((B(i,j))^2)+((R(i,j))^2)));
        eprx(i,j,9)=atan((R(i,j)/B(i,j)));
        
    end
end
 
for i=1:1:clno %first 2 randomly chosen points are the means mclu holds means of classes no of rows=no of classes, no of col=9
    xf=randi([1 row],1,1);
    yf=randi([1 col],1,1);
    clusters(1,1,i)=xf;
    clusters(1,2,i)=yf;
    for j=1:1:9
        
            mclu(i,j)=eprx(xf,yf,j);
        
    end
end
for t=1:1:10
m=zeros(1000,1);
 
 
temp=zeros(100,1);
 
    for i=1:1:row
        for j=1:1:col
            
            for k=1:1:clno
                s=0;
                for l=1:1:9
                   s1=(abs(eprx(i,j,l)-mclu(k,l)));
                    s=s+s1;
                    %s=1-exp(-s);
                end
                distance(k)=s;
            end
            [maxval,fclus]=min(distance);
            m(fclus)=m(fclus)+1;
            clusters(m(fclus),1,fclus)=i;
            clusters(m(fclus),2,fclus)=j;
            
       end
    end
    
     % after all points are done once calculate new means
    for os=1:1:clno
        for es=1:1:9
            sum=0;
            for is=1:1:m(os)
                xss=clusters(is,1,os);
                yss=clusters(is,2,os);
                sumo=eprx(xss,yss,es);
                sum=sum+sumo;
            end
            mclu(os,es)=sum/m(os);
        end
       red=randi(255,1);gr=randi(255,1);bl=randi(255,1);
        for color=1:1:m(os)%color all pixels in a cluster with the mean color
            xs1=clusters(color,1,os);
            ys1=clusters(color,2,os);
            
            y(xs1,ys1,1)=mclu(os,1);
            y(xs1,ys1,2)=mclu(os,2);
            y(xs1,ys1,3)=mclu(os,3);
        end
    end
end
                
  
 
 
%finding mean square error
error(:,:,:)=0;
error=(double(y) - double(x)) .^ 2;
suu=0;
for i=1:1:row
    for j=1:1:col
        for k=1:1:3
            su=(error(i,j,k));
            suu=su+suu;
        end
    end
end
 
ms=suu/(row*col*3);
 
%snr
sig=double(x).^2;
suusig=0;
for i=1:1:row
    for j=1:1:col
        for k=1:1:3
            susig=(sig(i,j,k));
            suusig=susig+suusig;
        end
    end
end
snr=10*log(suusig/suu);
