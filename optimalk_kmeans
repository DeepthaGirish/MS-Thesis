
clc;
clear all;
close all;
x1=imread('D:\University of Cincinnati\MS Thesis\Image dataset\1.jpg');
record(:,:)=0;
 
%y1=imread('C:\Users\Deeptha93\Desktop\night out230714\IMG_8331.JPG');
%x2=imread('C:\Users\Deeptha93\Documents\Visual Studio 2013\Projects\OpenCV-Test\OpenCV-Test\1.png');
 
 
 
    distance=zeros();
x= imresize(x1,.5);% this image is too big
 
big=size(x);
row=big(1);
col=big(2);
si=row*col;
y=zeros(row,col,3);
x=im2double(x);
eprx(:,:,:)=0;
 
ui=0;
clusters(:,:,:)=0;% each plane is one cluster and has point coordinates
mclu(:,:)=0;% holds mean vector for every cluster
R=x(:,:,1);G=x(:,:,2);B=x(:,:,3);% eprx is the new pixel rep as described in mam's paper
f=3;
for i=1:1:row
    for j=1:1:col
        eprx(i,j,1)=R(i,j);
        eprx(i,j,2)=G(i,j);
        eprx(i,j,3)=B(i,j);
       
    end
end
% 2d (pixel location) array of thikness k will be used to store k classes. A separate 1d
% array will be used to store the mean of each class
%k-means segmentation. k is not fixed.
 
clusters(1,1,1)=randi(row,1);% first pixel chosen
xf=clusters(1,1,1);
clusters(1,2,1)=randi(col,1);
yf=clusters(1,2,1);
clusters(1,1,2)=randi(row,1);% second pixel chosen
xs=clusters(1,1,2);
clusters(1,2,2)=randi(col,1);
ys=clusters(1,2,2);
clno=2;
 
countcl=ones(100,1);% number of elements in each cluster
 
for i=1:1:clno %first 2 randomly chosen points are the means mclu holds means of classes no of rows=no of classes, no of col=9
    for j=1:1:f
        if (i==1)
            mclu(i,j)=eprx(xf,yf,j);
        else
            mclu(i,j)=eprx(xs,ys,j);
        end
    end
end
 
 
for t=1:1:10
    m=zeros(1000,1);
    sumclus=zeros(1000,f);
    for i=1:1:row
        for j=1:1:col
            for k=1:1:clno
                %diatance=zeros(clno,1);
                s=0;
                for l=1:1:f
                   s1=((eprx(i,j,l)-mclu(k,l)).^2);
                    s=s+s1;
                    %s=1-exp(-s);
                end
                distance(k)=sqrt(s);
            end
            [minval,fclus]=min(distance);
            if(minval< (0.75*mean(distance))) % replace 0.75 with thresh
                
                m(fclus)=m(fclus)+1;
                clusters(m(fclus),1,fclus)=i;
                clusters(m(fclus),2,fclus)=j;
                for a1=1:1:f
                   sumclus(fclus,a1)=sumclus(fclus,a1)+eprx(i,j,a1);
                end
                
            else
                clno=clno+1;
                m(clno)=m(clno)+1;
                clusters(m(clno),1,clno)=i;
                clusters(m(clno),2,clno)=j;
                
                for q=1:1:f
                    mclu(clno,q)=eprx(i,j,q);
                end
            end
            end
       
    end
     for a2=1:1:clno
            for a3=1:1:f
               mclu(a2,a3)=sumclus(a2,a3)/m(a2);
            end
        end
 
end
clno,
 
% after all points are done once calculate new means
z=zeros(row,col);
for os=1:1:clno
    for es=1:1:f
        sum=0;
        for is=1:1:m(os)-1
            xss=clusters(is,1,os);
            yss=clusters(is,2,os);
            sumo=eprx(xss,yss,es);
            sum=sum+sumo;
        end
        mclu(os,es)=sum/m(os);
    end
    for color=1:1:m(os)%color all pixels in a cluster with the mean color
        xs1=clusters(color,1,os);
        ys1=clusters(color,2,os);
        z(xs1,ys1)=os;
        y(xs1,ys1,1)=mclu(os,1);
        y(xs1,ys1,2)=mclu(os,2);
        y(xs1,ys1,3)=mclu(os,3);
        
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
mse=suu/(si*3);
%%%% finding sil of all points and their cv
k=0;
for i=1:1:row
    for j=1:1:col
        k=k+1;
        for l=1:1:f
            
            data(k,l)=x(i,j,l);
            
           
        end
    end
end
 
k=0;ind=zeros();
for i=1:1:row
    for j=1:1:col
        k=k+1;
        ind(k)=z(i,j);
    end
end
s1=silhouette(data,ind);
cv1=mean(s1)/std(s1);
 
