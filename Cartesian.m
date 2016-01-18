paturn=1;
model;
bandPolygonNumber=100;         % バンド分割数
r=3.75;                 % バンド半径
P=0.45*10^(-7);         % ダイポールモーメント
eps=2.213545.*10.^(-4); % 誘電率

base.Pco=zeros(bandPolygonNumber,3);                            % バンドポリゴン座標
i=0:(2*pi)/bandPolygonNumber:2*pi*((bandPolygonNumber-1)/bandPolygonNumber); 

base.Pco(:,1)=r*cos(i);  % バンドポリゴン初期ｘ座標
base.Pco(:,3)=r*sin(i);  % バンドポリゴン初期ｚ座標
base.Pmo=-base.Pco/r*P;   % バンドポリゴン初期ダイポールモーメント
                        % base.mo(:,1) ｘ成分
                        % base.mo(:,3) ｚ成分

cel=-25:25;
map=zeros(numel(cel), numel(cel), 30);
map2=zeros(numel(cel), numel(cel), 30);
map3=zeros(numel(cel), numel(cel), 30);
Base=zeros(30,1);
Max=zeros(30,51);
MAX=zeros(1,30);

for time=1:30
    
    % 平行移動行列
    A=(modelpaturn(paturn).Oc(time,:))';
    A(4,1)=0;
    A1=cat(2,zeros(4,3),A)+eye(4);

    % ダミー行列の追加
    B=base.Pco';
    B1=vertcat(B,ones(1,bandPolygonNumber));

    % 平行移動
    Dposi=A1*B1;

    % ダミー行列の削除
    Dposi(4,:)=[];
    
    
    
    for tate=1:51
        for yoko=1:51
                       
            % ρの定義
            Eposi=[cel(yoko);
                   cel(tate);
                   6.75];
            Eposi=Eposi(:,ones(1,bandPolygonNumber));

            rho=Eposi-Dposi;
            map(tate,yoko,time)= trace((base.Pmo*rho)/(4*pi*eps*norm(rho)^3));
        end
    end
    
    Base(time)=map(26+20, 26-20, time);
    map2(:,:,time)=map(:,:,time)-Base(time);
    
    Max(time,:)=max(map2(:,:,time));
    MAX(time)=max(Max(time,:));
    map3(:,:,time)=map2(:,:,time)/MAX(time);
end

[X2,Y2]=meshgrid(cel,cel);
% figure;
for i=1:30
%     subplot(5,6,i);
    figure;
    surf(X2,Y2,map2(:,:,i));
    shading('flat');
%     colorbar;
    caxis([-MAX(i) MAX(i)]);
    xlim([-25 25]);
    ylim([-25 25]);
    set(gca,'XTick',[-25,-20,-15,-10,-5,0,5,10,15,20,25]);
    set(gca,'YTick',[-25,-20,-15,-10,-5,0,5,10,15,20,25]);
    view(0,90);
    
    name=strcat('C:\Users\m-kawano\Documents\参考\CTcolonoscopy\一時/',num2str(i));
    saveas(gcf, name, 'jpg')
end

%% avi出力
%{
obj=VideoWriter('C:\Users\m-kawano\Documents\参考\CTcolonoscopy\一時/modelname');
obj.FrameRate=1;
open(obj)
for i=1:30
    name=strcat('C:\Users\m-kawano\Documents\参考\CTcolonoscopy\一時/',num2str(i),'.jpg');
    image(imread(name));
    drawnow;
    writeVideo(obj,getframe);
end
close(obj);
%}