clear;
% --- �萔��` ---
Paturn=3;                                           % ���f���p�^�[���I��
model;                                              % ���f���Ăяo��
% Delta=11;                                           % �_�C�|�[���o���h��(11*10^-3cm)
% r=3.75;                                             % �_�C�|�[���o���h�W�����a(3.75cm)
Cel=-25:25;                                         % �̕\�ʔ͈�(50cm)
% Cel=-15:15;                                         % �~��
% Cel=-20:20;                                         % �~��
SampRate=modelpaturn(Paturn).SamplingRate;          % �o���h���_�̒��S���W��
DipoleBandRate=SampRate-1;                          % �_�C�|�[���o���h�̒��S���W��
% Polygon=floor(2*pi*r/Delta)*10^3;                   % �o���h�̒��_��(���Ȃ�)
Polygon=100;
Dipole=Polygon*2;                               % �_�C�|�[����
P=0.45*10^(-9);                                     % �_�C�|�[�����[�����g
% P=2.2*10^(-16);
% D=(0.45*10^-9)/(2*pi*1.25*0.9);
% D=10^-15*[0.0207 0.0296 0.0422 0.0590 0.0804 0.1044 0.1301 0.1573 0.1821 0.2040 0.2246 0.2317 0.2482 0.2420 0.2442 0.2292 0.2272 0.2455 0.4661];
Eps=2.213545*10^(-4);                             % �U�d��
% Eps=8.85418*10^(-12);
Theta=0:(2*pi)/Polygon:2*pi;  % �o���h���S�p

% --- �s���` ---
map =zeros(numel(Cel),numel(Cel), DipoleBandRate);              % �̕\�ʓd�ʕ��z(������)
map2=zeros(numel(Cel),numel(Cel), DipoleBandRate);              % �̕\�ʓd�ʕ��z(��d�ɍ���)
map3=zeros(numel(Cel),numel(Cel), DipoleBandRate);              % �̕\�ʓd�ʕ��z(���K��)
Base=zeros(DipoleBandRate,1);                                   % ��d�ɓd��
Max=zeros(DipoleBandRate,numel(Cel));                           % ����T���v���_�ł̍ő�d��(�s����)
MAX=zeros(1,DipoleBandRate);                                    % ����T���v���_�ł̍ő�d��
poly.Pos=zeros(3,Polygon+1,SampRate);                             % ���_�̏������W
poly.Mom=zeros(3,Polygon+1,SampRate);                             % ���_�̏������[�����g
poly.RotatP=zeros(3,Polygon+1,SampRate);                          % ���_�̕ϊ�����W
poly.TransP=zeros(4,Polygon+1,SampRate);                          % ���_�̕ϊ��ヂ�[�����g
dipole.Pos=zeros(3,Dipole,DipoleBandRate);                      % �_�C�|�[�����W
dipole.Mom=zeros(3,Dipole,DipoleBandRate);                      % �_�C�|�[�����[�����g
dipole.dS=zeros(Dipole,DipoleBandRate);                         % �����̈�
dis=zeros(1,Dipole);                                              % ��^3
dipole.dV=zeros(Dipole,numel(Cel),numel(Cel),DipoleBandRate);   % �P�̃_�C�|�[���ɂ��d��
rho=zeros(3,Dipole,numel(Cel),numel(Cel),DipoleBandRate);       % �_�C�|�[������d�ɂ܂ł̋���

% --- ���_���W ---
for srNum=1:SampRate;  % �o���h(���S���W)�ԍ�
    poly.Pos(1,:,srNum)=modelpaturn(Paturn).R(srNum)*cos(Theta);% ���_�̏��������W
    poly.Pos(3,:,srNum)=modelpaturn(Paturn).R(srNum)*sin(Theta);% ���_�̏��������W
    
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
    
    % ��]
    ZaxialRotation=[ cos(phi) sin(phi) 0;
                    -sin(phi) cos(phi) 0;
                        0       0      1];
    YaxialRotation=[ cos(theta) 0 sin(theta);
                        0       1      0    ;
                    -sin(theta) 0 cos(theta)];
    
    poly.RotatP(:,:,srNum)=YaxialRotation*ZaxialRotation*poly.Pos(:,:,srNum);
    
    
    % ���s�ړ��s��
    A=(modelpaturn(Paturn).BandCenter(srNum,:))';
    A(4,1)=0;
    A1=cat(2,zeros(4,3),A)+eye(4);

    % �_�~�[�s��̒ǉ�
    B=poly.RotatP(:,:,srNum);
    B1=vertcat(B,ones(1,Polygon+1));

    % ���s�ړ�
    poly.TransP(:,:,srNum)=A1*B1;    
end

% �_�~�[�s��̍폜
poly.TransP(4,:,:)=[];

    
% �_�C�|�[���ɂ��Ă��ꂼ��Z�o
for dbrNum=1:DipoleBandRate
    
    for pNum=1:Polygon
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
        dipole.Mom(:,2*pNum-1,dbrNum)=dipole.dS(2*pNum-1,dbrNum)*cross(nextupVec,upVec)/norm(cross(nextupVec,upVec));
        dipole.Mom(:,2*pNum,dbrNum)=dipole.dS(2*pNum,dbrNum)*cross(nextVec,nextupVec)/norm(cross(nextVec,nextupVec));        
    end
    
    D=P/sum(dipole.dS(:,dbrNum));
    dipole.Mom(:,:,dbrNum)=D*dipole.Mom(:,:,dbrNum);

end


for dbrNum=1:DipoleBandRate
    for yoko=1:numel(Cel)
        for tate=1:numel(Cel)
                       
            % �d�Ɉʒu
%             electrode=[Cel(yoko); Cel(tate); 6.75];
%             electrode=[Cel(yoko); 12; Cel(tate)]; % ����(�~��)
%             electrode=[Cel(yoko); 17; Cel(tate)]; % ����(�~��)
%             electrode=[Cel(yoko); Cel(tate); 15];   % z������
%             electrode=[Cel(yoko); 15; Cel(tate)];   % y������
            electrode=[15; Cel(tate); Cel(yoko)];   % x������
            elePos=electrode(:,ones(1,Dipole));

            % ��
            rho(:,:,tate,yoko,dbrNum)=elePos-dipole.Pos(:,:,dbrNum);
            for dipNum=1:Dipole
                dis(dipNum)=norm(rho(:,dipNum,tate,yoko,dbrNum));
            end
            
            % �_�C�|�[���ɂ��P�̓d�ɓd��
%             dipole.dV(:,tate,yoko,dbrNum)=diag(dipole.Mom(:,:,dbrNum)'*rho(:,:,tate,yoko,dbrNum))./(4*pi*Eps*dis(:).^3);
            map(tate,yoko,dbrNum)=sum(diag(-dipole.Mom(:,:,dbrNum)'*rho(:,:,tate,yoko,dbrNum))./(4*pi*Eps*dis'.^3));
            
            % �_�C�|�[���o���h�ɂ���̓d�ɓd��
%             map(tate,yoko,dbrNum)=(dipole.dV(:,tate,yoko,dbrNum)'*dipole.dS(:,dbrNum));
        end
    end
    
    Base(dbrNum)=map(26-20, 26+20, dbrNum);
%     Base(dbrNum)=map(16-10, 16+10, dbrNum); % �~��
%     Base(dbrNum)=map(21-10, 21+10, dbrNum); % �~��
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
%     subplot(3,4,i); % �~��
%     subplot(4,5,i); % �~��
%     figure;
%     surf(X2,Y2,map(:,:,i));
    surf(X2,Y2,map2(:,:,i));
    shading('flat');
%     colorbar;
%     caxis([-MAX(i) MAX(i)]);
%     caxis([-2*10^-10 4*10^-10]);    % �~��/����
%     caxis(10^-8*[-1 4]); % �~��1/����
%     caxis([0 8*10^-8]); % �~��2

    xlim([-25 25]); ylim([-25 25]);
%     xlim([-15 15]); ylim([-15 15]); % �~��
%     xlim([-20 20]); ylim([-20 20]); % �~��
    set(gca,'XTick',[-25,-20,-15,-10,-5,0,5,10,15,20,25]); set(gca,'YTick',[-25,-20,-15,-10,-5,0,5,10,15,20,25]);
    view(0,90);
    
%     name=strcat('C:\Users\m-kawano\Documents\�Q�l\CTcolonoscopy\�ꎞ/',num2str(i));
%     saveas(gcf, name, 'jpg')
end
%}
%% avi�o��
%{
obj=VideoWriter('C:\Users\m-kawano\Documents\�Q�l\CTcolonoscopy\�ꎞ/modelname');
obj.FrameRate=1;
open(obj)
for i=1:DipoleBandRate
    name=strcat('C:\Users\m-kawano\Documents\�Q�l\CTcolonoscopy\�ꎞ/',num2str(i),'.jpg');
    image(imread(name));
    drawnow;
    writeVideo(obj,getframe);
end
close(obj);
%}