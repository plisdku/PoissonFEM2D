%max(abs(xv_matrix(:,1) - VV.ParticleArray(1).xv))
figure(30)
imagesc((E_r(:,:,1) - VV.Er(:,:,1)./E_r(:,:,1)))

figure(4190)
imagesc(VV.Er(:,:,1)')

figure(5010)
imagesc(E_r(:,:,1)')

figure(500)
imagesc(VV.V(:,:,1)')

figure(510)
imagesc(U(:,:,1)')
%%
figure(566)
imagesc((VV.V(:,:,1)' - U(:,:,1)')./U(:,:,1)')

%%
figure(888)
imagesc((u_cartesian' - uu(:,:,1)')./u_cartesian')