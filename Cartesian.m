% --- 定数定義 ---
Paturn=1;                   % モデルパターン選択
model;                      % モデル呼び出し
Cel=-25:25;                 % 体表面範囲(50cm)
SampRate=30;                % バンド頂点の中心座標数
DipoleBandRate=SampRate-1;  % ダイポールバンドの中心座標数
Polygon=100;                % バンドの頂点数(幅なし)
Dipole=(Polygon-1)*2;       % ダイポール数
P=0.45*10^(-7);             % ダイポールモーメント
Eps=2.213545.*10.^(-4);     % 誘電率
Theta=0:(2*pi)/Polygon:2*pi*((Polygon-1)/Polygon);    % バンド中心角

% --- 行列定義 ---
map =zeros(numel(Cel),numel(Cel), DipoleBandRate);             % 体表面電位分布(無限遠)
map2=zeros(numel(Cel),numel(Cel), DipoleBandRate);             % 体表面電位分布(基準電極差分)
map3=zeros(numel(Cel),numel(Cel), DipoleBandRate);             % 体表面電位分布(正規化)
Base=zeros(DipoleBandRate,1);                                  % 基準電極電位
Max=zeros(DipoleBandRate,numel(Cel));                          % あるサンプル点での最大電位(行ごと)
MAX=zeros(1,DipoleBandRate);                                   % あるサンプル点での最大電位
poly.Pos=zeros(3,Polygon,SampRate);                         % 頂点の初期座標
poly.Mom=zeros(3,Polygon,SampRate);                         % 頂点の初期モーメント
poly.TransP=zeros(4,Polygon,SampRate);                      % 頂点の変換後座標
poly.TransM=zeros(3,Polygon,SampRate);                      % 頂点の変換後モーメント
dipole.Cen=zeros(3,Dipole,DipoleBandRate);               % ダイポールバンドの中心座標
dipole.Pos=zeros(3,Dipole,DipoleBandRate);               % ダイポール座標
dipole.Mom=zeros(3,Dipole,DipoleBandRate);               % ダイポールモーメント
dipole.dS=zeros(Dipole,DipoleBandRate);                  % 微小領域
dipole.dV=zeros(Dipole,numel(Cel),numel(Cel),DipoleBandRate);           % １つのダイポールによる電位
rho=zeros(3,Dipole,numel(Cel),numel(Cel),DipoleBandRate);        % ダイポールから電極までの距離

% --- 頂点座標 ---
for srNum=1:SampRate;  % バンド(中心座標)番号
    poly.Pos(1,:,srNum)=modelpaturn(Paturn).R(srNum)*cos(Theta);                % 頂点の初期ｘ座標
    poly.Pos(3,:,srNum)=modelpaturn(Paturn).R(srNum)*sin(Theta);                % 頂点の初期ｚ座標
    % 平行移動行列
    A=(modelpaturn(Paturn).BandCenter(srNum,:))';
    A(4,1)=0;
    A1=cat(2,zeros(4,3),A)+eye(4);

    % ダミー行列の追加
    B=poly.Pos(:,:,srNum);
    B1=vertcat(B,ones(1,Polygon));

    % 平行移動
    poly.TransP(:,:,srNum)=A1*B1;

    % 方向変換
    
end

% ダミー行列の削除
poly.TransP(4,:,:)=[];
    

for dbrNum=1:DipoleBandRate
    for pNum=1:Polygon-1
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

    end
    
    dC=(modelpaturn(Paturn).BandCenter(dbrNum,:)+modelpaturn(Paturn).BandCenter(dbrNum+1,:))/2;
    dipole.Cen(:,:,dbrNum)=repmat(dC',1,Dipole);
end
    
dipole.Mom=dipole.Cen-dipole.Pos;

for dbrNum=1:DipoleBandRate
    for yoko=1:numel(Cel)
        for tate=1:numel(Cel)
                       
            % ρの定義
            C=[Cel(yoko);
               Cel(tate);
               6.75];
            elePos=C(:,ones(1,Dipole));

            rho(:,:,tate,yoko,dbrNum)=elePos-dipole.Pos(:,:,dbrNum);
            dipole.dV(:,tate,yoko,dbrNum)=diag((dipole.Mom(:,:,dbrNum)'*rho(:,:,tate,yoko,dbrNum))/(4*pi*Eps*norm(rho(:,:,tate,yoko,dbrNum))^3));
            map(tate,yoko,dbrNum)=dipole.dV(:,tate,yoko,dbrNum)'*dipole.dS(:,dbrNum);
        end
    end
    
    Base(dbrNum)=map(26+20, 26-20, dbrNum);
    map2(:,:,dbrNum)=map(:,:,dbrNum)-Base(dbrNum);
    
    Max(dbrNum,:)=max(map2(:,:,dbrNum));
    MAX(dbrNum)=max(Max(dbrNum,:));
    map3(:,:,dbrNum)=map2(:,:,dbrNum)/MAX(dbrNum);
end

[X2,Y2]=meshgrid(Cel,Cel);
% figure;
for i=1:DipoleBandRate
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
    
%     name=strcat('C:\Users\m-kawano\Documents\参考\CTcolonoscopy\一時/',num2str(i));
%     saveas(gcf, name, 'jpg')
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