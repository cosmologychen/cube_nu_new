use omp_lib
use parameters
implicit none
real,allocatable :: Gk1(:,:,:),Gk2(:,:,:),Gk3(:,:,:) 

call geometry
call omp_set_num_threads(ncore)
if (head) then
  print*,'kernel_initialization on',int(ncore,kind=2),' cores'
  print*,'  nc,ngt,nft =',nc,ngt,nft
endif

allocate(Gk1(nc/2+1,nc,nc),Gk2(ngt/2+1,ngt,ngt),Gk3(nft/2+1,nft,nft))
call Green_3D(Gk1,nc,nc/2+1,nc,nc,apm1c,0.,real(ratio_cs))
call Green_3D(Gk2,ngt,ngt/2+1,ngt,ngt,apm2,apm1,1.)
call Green_3D(Gk3,nft,nft/2+1,nft,nft,apm3f,apm2f,1./ratio_sf)

call system('mkdir -p '//opath//'image'//image2str(image))
open(11,file=output_dir()//'Gk1'//output_suffix(),access='stream')
write(11) Gk1
close(11)
if (head) then
  open(11,file=output_dir()//'Gk2.bin',access='stream')
  write(11) Gk2
  close(11)
  open(11,file=output_dir()//'Gk3.bin',access='stream')
  write(11) Gk3
  close(11)
endif

deallocate(Gk1,Gk2,Gk3)

end