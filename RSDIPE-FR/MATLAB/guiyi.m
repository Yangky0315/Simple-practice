function [ output ] = guiyi( input )
input=double(input);
[x,y,z]=size(input);

for k=1:z
    maxi=max(max(input(:,:,k)));
    mini=min(min(input(:,:,k)));
    for i=1:x
        for j=1:y
            %πÈ“ªªØ
            output(i,j,k)=(input(i,j,k)-mini)/(maxi-mini);
        end
    end
end
end

