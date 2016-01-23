% --- �萔��` ---
Paturn=1;                   % ���f���p�^�[���I��
model;                      % ���f���Ăяo��
Cel=-25:25;                 % �̕\�ʔ͈�(50cm)
SampRate=30;                % �o���h���_�̒��S���W��
DipoleBandRate=SampRate-1;  % �_�C�|�[���o���h�̒��S���W��
Polygon=100;                % �o���h�̒��_��(���Ȃ�)
Dipole=(Polygon-1)*2;       % �_�C�|�[����
P=0.45*10^(-7);             % �_�C�|�[�����[�����g
Eps=2.213545.*10.^(-4);     % �U�d��
Theta=0:(2*pi)/Polygon:2*pi*((Polygon-1)/Polygon);    % �o���h���S�p

% --- �s���` ---
map =zeros(numel(Cel),numel(Cel), DipoleBandRate);             % �̕\�ʓd�ʕ��z(������)
map2=zeros(numel(Cel),numel(Cel), DipoleBandRate);             % �̕\�ʓd�ʕ��z(��d�ɍ���)
map3=zeros(numel(Cel),numel(Cel), DipoleBandRate);             % �̕\�ʓd�ʕ��z(���K��)
Base=zeros(DipoleBandRate,1);                                  % ��d�ɓd��
Max=zeros(DipoleBandRate,numel(Cel));                          % ����T���v���_�ł̍ő�d��(�s����)
MAX=zeros(1,DipoleBandRate);                                   % ����T���v���_�ł̍ő�d��
poly.Pos=zeros(3,Polygon,SampRate);                         % ���_�̏������W
poly.Mom=zeros(3,Polygon,SampRate);                         % ���_�̏������[�����g
poly.TransP=zeros(4,Polygon,SampRate);                      % ���_�̕ϊ�����W
poly.TransM=zeros(3,Polygon,SampRate);                      % ���_�̕ϊ��ヂ�[�����g
dipole.Cen=zeros(3,Dipole,DipoleBandRate);               % �_�C�|�[���o���h�̒��S���W
dipole.Pos=zeros(3,Dipole,DipoleBandRate);               % �_�C�|�[�����W
dipole.Mom=zeros(3,Dipole,DipoleBandRate);               % �_�C�|�[�����[�����g
dipole.dS=zeros(Dipole,DipoleBandRate);                  % �����̈�
dipole.dV=zeros(Dipole,numel(Cel),numel(Cel),DipoleBandRate);           % �P�̃_�C�|�[���ɂ��d��
rho=zeros(3,Dipole,numel(Cel),numel(Cel),DipoleBandRate);        % �_�C�|�[������d�ɂ܂ł̋���

% --- ���_���W ---
for srNum=1:SampRate;  % �o���h(���S���W)�ԍ�
    poly.Pos(1,:,srNum)=modelpaturn(Paturn).R(srNum)*cos(Theta);                % ���_�̏��������W
    poly.Pos(3,:,srNum)=modelpaturn(Paturn).R(srNum)*sin(Theta);                % ���_�̏��������W
    % ���s�ړ��s��
    A=(modelpaturn(Paturn).BandCenter(srNum,:))';
    A(4,1)=0;
    A1=cat(2,zeros(4,3),A)+eye(4);

    % �_�~�[�s��̒ǉ�
    B=poly.Pos(:,:,srNum);
    B1=vertcat(B,ones(1,Polygon));

    % ���s�ړ�
    poly.TransP(:,:,srNum)=A1*B1;

    % �����ϊ�
    
end

% �_�~�[�s��̍폜
poly.TransP(4,:,:)=[];
    

for dbrNum=1:DipoleBandRate
    for pNum=1:Polygon-1
        % �_�C�|�[�����W
        dipole.Pos(:,2*pNum-1,dbrNum)=(poly.TransP(:,pNum,dbrNum)+poly.TransP(:,pNum,dbrNum+1)+poly.TransP(:,pNum+1,dbrNum+1))/3;
        dipole.Pos(:,2*pNum,dbrNum)=(poly.TransP(:,pNum,dbrNum)+poly.TransP(:,pNum+1,dbrNum)+poly.TransP(:,pNum+1,dbrNum+1))/3;
        
        % �����̈�
        upVec=poly.TransP(:,pNum,dbrNum+1)-poly.TransP(:,pNum,dbrNum);
        nextVec=poly.TransP(:,pNum+1,dbrNum)-poly.TransP(:,pNum,dbrNum);
        nextupVec=poly.TransP(:,pNum+1,dbrNum+1)-poly.TransP(:,pNum,dbrNum);
        dipole.dS(2*pNum-1,dbrNum)=norm(cross(upVec,nextupVec))/2;
        dipole.dS(2*pNum,dbrNum)=norm(cross(nextVec,nextupVec))/2;
        
        % �_�C�|�[�����[�����g

    end
    
    dC=(modelpaturn(Paturn).BandCenter(dbrNum,:)+modelpaturn(Paturn).BandCenter(dbrNum+1,:))/2;
    dipole.Cen(:,:,dbrNum)=repmat(dC',1,Dipole);
end
    
dipole.Mom=dipole.Cen-dipole.Pos;

for dbrNum=1:DipoleBandRate
    for yoko=1:numel(Cel)
        for tate=1:numel(Cel)
                       
            % �ς̒�`
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
    
%     name=strcat('C:\Users\m-kawano\Documents\�Q�l\CTcolonoscopy\�ꎞ/',num2str(i));
%     saveas(gcf, name, 'jpg')
end

%% avi�o��
%{
obj=VideoWriter('C:\Users\m-kawano\Documents\�Q�l\CTcolonoscopy\�ꎞ/modelname');
obj.FrameRate=1;
open(obj)
for i=1:30
    name=strcat('C:\Users\m-kawano\Documents\�Q�l\CTcolonoscopy\�ꎞ/',num2str(i),'.jpg');
    image(imread(name));
    drawnow;
    writeVideo(obj,getframe);
end
close(obj);
%}