function HamMPO = MPO_XXZ_full(L,J,Delta)
%MPO_XXZ_FULL Summary of this function goes here
%   Detailed explanation goes here
d = 2;
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0 ];
sigma_z = [1 0; 0 -1];
if L<3
  error('L should be at least 3');
end

HamMPO = cell(1,L);
HamMPO{1} = zeros(d,d,1,5);
HamMPO{1}(:,:,1,1) = zeros(d);
HamMPO{1}(:,:,1,2) = eye(d);
HamMPO{1}(:,:,1,3) = sigma_x;
HamMPO{1}(:,:,1,4) = 1i*sigma_y;   %i is introduced to make the MPO tensor real
HamMPO{1}(:,:,1,5) = sigma_z;
for pos = 2:(L-1)
  HamMPO{pos} = zeros(d,d,5,5);
  HamMPO{pos}(:,:,1,1) = eye(d);
  HamMPO{pos}(:,:,2,1) = zeros(d);
  HamMPO{pos}(:,:,3,1) = J*sigma_x;
  HamMPO{pos}(:,:,4,1) = -1i*J*sigma_y;  %i is introduced to make the MPO tensor real
  HamMPO{pos}(:,:,5,1) = J*Delta*sigma_z;
  HamMPO{pos}(:,:,2,2) = eye(d);
  HamMPO{pos}(:,:,2,3) = sigma_x;
  HamMPO{pos}(:,:,2,4) = 1i*sigma_y;  %i is introduced to make the MPO tensor real
  HamMPO{pos}(:,:,2,5) = sigma_z;
end
HamMPO{L} = zeros(d,d,5,1);
HamMPO{L}(:,:,1,1) = eye(d);
HamMPO{L}(:,:,2,1) = zeros(d);
HamMPO{L}(:,:,3,1) = J*sigma_x;
HamMPO{L}(:,:,4,1) = -1i*J*sigma_y;  %i is introduced to make the MPO tensor real
HamMPO{L}(:,:,5,1) = J*Delta*sigma_z;

end

