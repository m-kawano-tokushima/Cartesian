clear;
% --- 定数定義 ---
Paturn=3;                                           % モデルパターン選択
model;                                              % モデル呼び出し
% Delta=11;                                           % ダイポールバンド幅(11*10^-3cm)
% r=3.75;                                             % ダイポールバンド標準半径(3.75cm)
Cel=-25:25;                                         % 体表面範囲(50cm)
% Cel=-15:15;                                         % 円柱
% Cel=-20:20;                                         % 円錐
SampRate=modelpaturn(Paturn).SamplingRate;          % バンド頂点の中心座標数
DipoleBandRate=SampRate-1;                          % ダイポールバンドの中心座標数
% Polygon=floor(2*pi*r/Delta)*10^3;                   % バンドの頂点数(幅なし)
Polygon=100;
Dipole=Polygon*2;                               % ダイポール数
P=0.45*10^(-9);                                     % ダイポールモーメント
% P=2.2*10^(-16);
% D=(0.45*10^-9)/(2*pi*1.25*0.9);
% D=10^-15*[0.0207 0.0296 0.0422 0.0590 0.0804 0.1044 0.1301 0.1573 0.1821 0.2040 0.2246 0.2317 0.2482 0.2420 0.2442 0.2292 0.2272 0.2455 0.4661];
Eps=2.213545*10^(-4);                             % 誘電率
% Eps=8.85418*10^(-12);
Theta=0:(2*pi)/Polygon:2*pi;  % バンド中心角

% --- 行列定義 ---
map =zeros(numel(Cel),numel(Cel), DipoleBandRate);              % 体表面電位分布(無限遠)
map2=zeros(numel(Cel),numel(Cel), DipoleBandRate);              % 体表面電位分布(基準電極差分)
map3=zeros(numel(Cel),numel(Cel), DipoleBandRate);              % 体表面電位分布(正規化)
Base=zeros(DipoleBandRate,1);                                   % 基準電極電位
Max=zeros(DipoleBandRate,numel(Cel));                           % あるサンプル点での最大電位(行ごと)
MAX=zeros(1,DipoleBandRate);                                    % あるサンプル点での最大電位
poly.Pos=zeros(3,Polygon+1,SampRate);                             % 頂点の初期座標
poly.Mom=zeros(3,Polygon+1,SampRate);                             % 頂点の初期モーメント
poly.RotatP=zeros(3,Polygon+1,SampRate);                          % 頂点の変換後座標
poly.TransP=zeros(4,Polygon+1,SampRate);                          % 頂点の変換後モーメント
dipole.Pos=zeros(3,Dipole,DipoleBandRate);                      % ダイポール座標
dipole.Mom=zeros(3,Dipole,DipoleBandRate);                      % ダイポールモーメント
dipole.dS=zeros(Dipole,DipoleBandRate);                         % 微小領域
dis=zeros(1,Dipole);                                              % ρ^3
dipole.dV=zeros(Dipole,numel(Cel),numel(Cel),DipoleBandRate);   % １つのダイポールによる電位
rho=zeros(3,Dipole,numel(Cel),numel(Cel),DipoleBandRate);       % ダイポールから電極までの距離

% --- 頂点座標 ---
for srNum=1:SampRate;  % バンド(中心座標)番号
    poly.Pos(1,:,srNum)=modelpaturn(Paturn).R(srNum)*cos(Theta);% 頂点の初期ｘ座標
    poly.Pos(3,:,srNum)=modelpaturn(Paturn).R(srNum)*sin(Theta);% 頂点の初期ｚ座標
    
    A=modelpaturn(Paturn).Moment(srNum,:)';
    
    if isnan(A(3)/A(1))
        theta=0;
    elseif A(3)==0
            if A(1)>0
                theta=0;
            else
                theta=pi;
            end
    else
        theta=-atan(A(3)/A(1));
    end
    phi=acos(A(2)/norm(A));
    
    % 回転
    ZaxialRotation=[ cos(phi) sin(phi) 0;
                    -sin(phi) cos(phi) 0;
                        0       0      1];
    YaxialRotation=[ cos(theta) 0 sin(theta);
                        0       1      0    ;
                    -sin(theta) 0 cos(theta)];
    
    poly.RotatP(:,:,srNum)=YaxialRotation*ZaxialRotation*poly.Pos(:,:,srNum);
    
    
    % 平行移動行列
    A=(modelpaturn(Paturn).BandCenter(srNum,:))';
    A(4,1)=0;
    A1=cat(2,zeros(4,3),A)+eye(4);

    % ダミー行列の追加
    B=poly.RotatP(:,:,srNum);
    B1=vertcat(B,ones(1,Polygon+1));

    % 平行移動
    poly.TransP(:,:,srNum)=A1*B1;    
end

% ダミー行列の削除
poly.TransP(4,:,:)=[];

    
% ダイポールについてそれぞれ算出
for dbrNum=1:DipoleBandRate
    
    for pNum=1:Polygon
        % ダイポール座標
        dipole.Pos(:,2*pNum-1,dbrNum)=(poly.TransP(:,pNum,dbrNum)+poly.TransP(:,pNum,dbrNum+1)+poly.TransP(:,pNum+1,dbrNum+1))/3;
        dipole.Pos(:,2*pNum,dbrNum)=(poly.TransP(:,pNum,dbrNum)+poly.TransP(:,pNum+1,dbrNum)+poly.TransP(:,pNum+1,dbrNum+1))/3;       

        % 微小領域
        upVec=poly.TransP(:,pNum,dbrNum+1)-poly.TransP(:,pNum,dbrNum);
        nextVec=poly.TransP(:,pNum+1,dbrNum)-poly.TransP(:,pNum,dbrNum);
        nextupVec=poly.TransP(:,pNum+1,dbrNum+1)-poly.TransP(:,pNum,dbrNum);
        dipole.dS(2*pNum-1,dbrNum)=norm(cross(upVec,nextupVec))/2;
        dipole.dS(2*pNum,dbrNum)=norm(cross(nextVec,nextupVec))/2;
                
        % ダイポールモーメント    
        dipole.Mom(:,2*pNum-1,dbrNum)=dipole.dS(2*pNum-1,dbrNum)*cross(nextupVec,upVec)/norm(cross(nextupVec,upVec));
        dipole.Mom(:,2*pNum,dbrNum)=dipole.dS(2*pNum,dbrNum)*cross(nextVec,nextupVec)/norm(cross(nextVec,nextupVec));        
    end
    
    D=P/sum(dipole.dS(:,dbrNum));
    dipole.Mom(:,:,dbrNum)=D*dipole.Mom(:,:,dbrNum);

end


for dbrNum=1:DipoleBandRate
    for yoko=1:numel(Cel)
        for tate=1:numel(Cel)
                       
            % 電極位置
%             electrode=[Cel(yoko); Cel(tate); 6.75];
%             electrode=[Cel(yoko); 12; Cel(tate)]; % 俯瞰(円柱)
%             electrode=[Cel(yoko); 17; Cel(tate)]; % 俯瞰(円錐)
%             electrode=[Cel(yoko); Cel(tate); 15];   % z軸平面
%             electrode=[Cel(yoko); 15; Cel(tate)];   % y軸平面
            electrode=[15; Cel(tate); Cel(yoko)];   % x軸平面
            elePos=electrode(:,ones(1,Dipole));

            % ρ
            rho(:,:,tate,yoko,dbrNum)=elePos-dipole.Pos(:,:,dbrNum);
            for dipNum=1:Dipole
                dis(dipNum)=norm(rho(:,dipNum,tate,yoko,dbrNum));
            end
            
            % ダイポールによる１つの電極電位
%             dipole.dV(:,tate,yoko,dbrNum)=diag(dipole.Mom(:,:,dbrNum)'*rho(:,:,tate,yoko,dbrNum))./(4*pi*Eps*dis(:).^3);
            map(tate,yoko,dbrNum)=sum(diag(-dipole.Mom(:,:,dbrNum)'*rho(:,:,tate,yoko,dbrNum))./(4*pi*Eps*dis'.^3));
            
            % ダイポールバンドによる一つの電極電位
%             map(tate,yoko,dbrNum)=(dipole.dV(:,tate,yoko,dbrNum)'*dipole.dS(:,dbrNum));
        end
    end
    
    Base(dbrNum)=map(26-20, 26+20, dbrNum);
%     Base(dbrNum)=map(16-10, 16+10, dbrNum); % 円柱
%     Base(dbrNum)=map(21-10, 21+10, dbrNum); % 円錐
    map2(:,:,dbrNum)=map(:,:,dbrNum)-Base(dbrNum);
    
    Max(dbrNum,:)=max(map2(:,:,dbrNum));
    MAX(dbrNum)=max(Max(dbrNum,:));
    map3(:,:,dbrNum)=map2(:,:,dbrNum)/MAX(dbrNum);
end
% %{
[X2,Y2]=meshgrid(Cel,Cel);
% figure;
for i=1:DipoleBandRate
    subplot(10,10,i);
%     subplot(3,4,i); % 円柱
%     subplot(4,5,i); % 円錐
%     figure;
%     surf(X2,Y2,map(:,:,i));
    surf(X2,Y2,map2(:,:,i));
    shading('flat');
%     colorbar;
%     caxis([-MAX(i) MAX(i)]);
%     caxis([-2*10^-10 4*10^-10]);    % 円柱/俯瞰
%     caxis(10^-8*[-1 4]); % 円錐1/俯瞰
%     caxis([0 8*10^-8]); % 円錐2

    xlim([-25 25]); ylim([-25 25]);
%     xlim([-15 15]); ylim([-15 15]); % 円柱
%     xlim([-20 20]); ylim([-20 20]); % 円錐
    set(gca,'XTick',[-25,-20,-15,-10,-5,0,5,10,15,20,25]); set(gca,'YTick',[-25,-20,-15,-10,-5,0,5,10,15,20,25]);
    view(0,90);
    
%     name=strcat('C:\Users\m-kawano\Documents\参考\CTcolonoscopy\一時/',num2str(i));
%     saveas(gcf, name, 'jpg')
end
%}
%% avi出力
%{
obj=VideoWriter('C:\Users\m-kawano\Documents\参考\CTcolonoscopy\一時/modelname');
obj.FrameRate=1;
open(obj)
for i=1:DipoleBandRate
    name=strcat('C:\Users\m-kawano\Documents\参考\CTcolonoscopy\一時/',num2str(i),'.jpg');
    image(imread(name));
    drawnow;
    writeVideo(obj,getframe);
end
close(obj);
%}