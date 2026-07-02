program ft8d

! Decode FT8 data read from *.c2 files or 12 kHz 16-bit PCM mono *.wav files.

  use iso_fortran_env, only: int8, int16, int32, error_unit
  include 'ft8_params.f90'
  character infile*80,msg37*37,date*6,time*4
  character msgcall*13,msggrid*4
  character*37 seenmessages(MAXCAND)
  real s(NFFT1,NHSYM)
  real sbase(NFFT1)
  real candidate(3,MAXCAND)
  real*8 dialfreq
  complex dd(NMAX,4)
  logical newdat,lft8apon,lsubtract,ldupe,is_wav
  real*8 seenxdt(MAXCAND), seenf1(MAXCAND)
  integer nseen
  integer apsym(58)
  integer wav_ierr,nparts,nwavcfg,iwcfg
  integer env_status,nsync_calls,nft8b_calls
  character(len=8) prof_env
  logical profile_on
  real t0,t1,t_read,t_sync,t_ft8b

  real :: elapsed, elapsed_start

  call cpu_time(elapsed_start)

  apsym=0
  apsym(1)=99
  apsym(30)=99

  prof_env=' '
  call get_environment_variable('FT8D_PROFILE',prof_env,status=env_status)
  profile_on = (env_status.eq.0 .and. len_trim(prof_env).gt.0 .and. prof_env(1:1).ne.'0')
  t_read=0.0
  t_sync=0.0
  t_ft8b=0.0
  nsync_calls=0
  nft8b_calls=0
  nseen=0
  seenmessages='                                     '
  seenxdt=0.0d0
  seenf1=0.0d0

  nargs=iargc()
  if(nargs.ne.1) then
    print*,'Usage: ft8d <file.c2|file.wav>'
    stop 1
  endif

  call getarg(1,infile)
  is_wav = index(infile,'.wav').gt.0
  
  dd = cmplx(0.0,0.0)
  if(is_wav) then
    if(profile_on) call cpu_time(t0)
    call read_wav_12k_mono(infile,dd(1,1),dialfreq,wav_ierr,.true.)
    if(wav_ierr.ne.0) stop 1
    call read_wav_12k_mono(infile,dd(1,2),dialfreq,wav_ierr,.false.)
    if(wav_ierr.ne.0) stop 1
    if(profile_on) then
      call cpu_time(t1)
      t_read=t_read + (t1-t0)
    endif
    date='000000'
    time='0000'
    nfa=-1800
    nfb=+1600
    nparts=2
  else
    nfa=-1600
    nfb=+1600
    nparts=4
    if(profile_on) call cpu_time(t0)
    open(10,file=infile,status='old',access='stream',form='unformatted')
    read(10) dialfreq,dd
    close(10)
    if(profile_on) then
      call cpu_time(t1)
      t_read=t_read + (t1-t0)
    endif
    j2=index(infile,'.c2')
    if(j2.ge.12) then
      date=infile(j2-11:j2-6)
      time=infile(j2-4:j2-1)
    else
      date='000000'
      time='0000'
    endif
  endif

  nfqso=0

  write(*,*) 'nparts ',nparts, ' npass=',npass, ' dialfreq=',dialfreq,' Hz, nfa=',nfa,' nfb=',nfb,' Hz'

  do ipart=1,nparts !4 for .c2 and 2 for .wav
    ndecodes=0
    if(is_wav) then
      nwavcfg=2
    else
      nwavcfg=1
    endif

    do iwcfg=1,nwavcfg  ! 1 for .c2 and 2 for .wav
      nQSOProgress=0
      n2=0
      ncontest=0
      lft8apon=.false.

      if(is_wav) then
        if(iwcfg.eq.1) then
          ndepth=3
        else
          ndepth=2
        endif
        syncmin=-3.0
      else
        ndepth=1
        syncmin=1.5
      endif

      if(ndepth.eq.1) npass=1
      if(ndepth.eq.2) npass=3
      if(ndepth.ge.3) npass=5

      do ipass=1,npass
        newdat=.true.
        if(.not.is_wav) syncmin=1.5

        ! Determine whether to subtract the previous pass's decodes 
        ! from the current pass's decodes.
        if(ipass.eq.1) then
          lsubtract=.true.
          if(ndepth.eq.1) lsubtract=.false.
        elseif(ipass.eq.2) then
          n2=ndecodes
          if(ndecodes.eq.0) cycle
          lsubtract=.true.
        elseif(ipass.eq.3) then
          if((ndecodes-n2).eq.0) cycle
          lsubtract=.false.
        endif

        if(profile_on) call cpu_time(t0)
        ! Call sync8 to find candidate decodes in the current pass.
        call sync8(dd(1:NMAX,ipart),nfa+2000,nfb+2000,syncmin, &
          nfqso+2000,s,candidate,ncand,sbase)
        if(profile_on) then
          call cpu_time(t1)
          t_sync=t_sync + (t1-t0)
          nsync_calls=nsync_calls+1
        endif

        do icand=1,ncand
          sync=candidate(3,icand)
          f1=candidate(1,icand)
          xdt=candidate(2,icand)
          xbase=10.0**(0.1*(sbase(nint(f1/3.125))-40.0))

          if(profile_on) call cpu_time(t0)
          ! Call ft8b to decode every candidate.
          ! write(*,*) 'Decoding candidate ',icand,' of ',ncand,' (pass ',ipass,')'
          ! write(*,*) '  NMAX=',NMAX,' f1=',f1,' xdt=',xdt,' sec, sync=',sync,' dB','ndepth=',ndepth,' lsubtract=',lsubtract
          call ft8b(dd(1:NMAX,ipart),newdat,nQSOProgress,nfqso+2000, &
            nftx,ndepth,lft8apon,lapcqonly,napwid,lsubtract,nagain, &
            ncontest,iaptype,f1,xdt,xbase,apsym,nharderrors,dmin, &
            nbadcrc,iappass,msg37,msgcall,msggrid,xsnr)
          if(profile_on) then
            call cpu_time(t1)
            t_ft8b=t_ft8b + (t1-t0)
            nft8b_calls=nft8b_calls+1
          endif

          nsnr=nint(xsnr)
          xdt=xdt-0.5
          hd=nharderrors+dmin

          if(nbadcrc.eq.0) then
            if(.not.is_wav) then
              if(msgcall(1:1).eq.' '.or.msgcall(1:1).eq.'<') cycle
            endif

            ldupe=.false.
            do id=1,nseen
                  if(msg37.eq.seenmessages(id) .and. &
                    abs(xdt-seenxdt(id)).le.0.5d0 .and. &
                    abs(f1-seenf1(id)).le.5.0d0) ldupe=.true.
            enddo
            if(ldupe) cycle

            ndecodes=ndecodes+1
            nseen=nseen+1
            seenmessages(nseen)=msg37
            seenxdt(nseen)=xdt
            seenf1(nseen)=f1

            call cpu_time(elapsed)
            
            if(is_wav) then
              write(*,2004) elapsed, date,nsnr,xdt,nint(f1-2000+dialfreq),msg37
2004        format(f10.6,1x,a6,1x,i4,f5.1,i5,1x,'~  ',a37)
            else
              write(*,1004) elapsed, date,time,15*(ipart-1),min(sync,999.0),nint(xsnr), &
                  xdt,nint(f1-2000+dialfreq),msgcall,msggrid
1004        format(f10.6,1x,a6,1x,a4,i2.2,f6.1,i4,f6.2,i9,1x,a13,1x,a4)
            endif
          endif
        enddo
      enddo
    enddo
  enddo       ! ipart=1,nparts

  call cpu_time(elapsed)
  write(*,*) 'Elapsed time (seconds): ',elapsed-elapsed_start

  if(profile_on) then
    write(error_unit,'(a,f8.4)') 'PROFILE read_s=',t_read
    write(error_unit,'(a,f8.4)') 'PROFILE sync_s=',t_sync
    write(error_unit,'(a,f8.4)') 'PROFILE ft8b_s=',t_ft8b
    write(error_unit,'(a,i0)') 'PROFILE sync_calls=',nsync_calls
    write(error_unit,'(a,i0)') 'PROFILE ft8b_calls=',nft8b_calls
  endif

contains

  subroutine read_wav_12k_mono(fname,ddout,dialfreq,ierr,use_fir)
    use iso_fortran_env, only: int16, int32
    character(len=*), intent(in) :: fname
    complex, intent(out) :: ddout(NMAX)
    real*8, intent(out) :: dialfreq
    integer, intent(out) :: ierr
    logical, intent(in) :: use_fir

    integer :: u, ios, pos0, i, k, idx, nout, it
    character(len=4) :: riff, wave, chunk
    integer(int32) :: riff_size, chunk_size, sample_rate, byte_rate, data_size
    integer(int16) :: audio_fmt, num_channels, block_align, bits_per_sample
    integer :: nframes
    integer(int16), allocatable :: pcm(:)
    real :: xr
    real*8 :: pi, dphi, phi
    real*8 :: fs, fc, xw, hk, win, hsum, sumr, sumi
    integer, parameter :: ntap=25
    integer, parameter :: m=(ntap-1)/2
    real*8 :: h(ntap)

    ierr = 0
    ddout = cmplx(0.0,0.0)
    dialfreq = 0.0d0

    open(newunit=u,file=fname,status='old',access='stream',form='unformatted',action='read',iostat=ios)
    if(ios.ne.0) then
      print*,'Error opening WAV file: ',trim(fname)
      ierr = 1
      return
    endif

    read(u,iostat=ios) riff
    read(u,iostat=ios) riff_size
    read(u,iostat=ios) wave
    if(ios.ne.0 .or. riff.ne.'RIFF' .or. wave.ne.'WAVE') then
      print*,'Unsupported WAV header in ',trim(fname)
      close(u)
      ierr = 2
      return
    endif

    sample_rate = 0
    data_size = 0

    do
      read(u,iostat=ios) chunk
      if(ios.ne.0) exit
      read(u,iostat=ios) chunk_size
      if(ios.ne.0) exit

      if(chunk.eq.'fmt ') then
        read(u) audio_fmt
        read(u) num_channels
        read(u) sample_rate
        read(u) byte_rate
        read(u) block_align
        read(u) bits_per_sample
        if(chunk_size.gt.16) then
          inquire(u,pos=pos0)
          read(u,pos=pos0 + (chunk_size-16))
        endif
      elseif(chunk.eq.'data') then
        data_size = chunk_size
        exit
      else
        inquire(u,pos=pos0)
        read(u,pos=pos0 + chunk_size)
      endif
    enddo

    if(data_size.le.0) then
      print*,'No data chunk in WAV: ',trim(fname)
      close(u)
      ierr = 3
      return
    endif

    if(audio_fmt.ne.1 .or. num_channels.ne.1 .or. bits_per_sample.ne.16 .or. sample_rate.ne.12000) then
      print*,'Unsupported WAV format (need PCM mono 16-bit 12kHz): ',trim(fname)
      close(u)
      ierr = 4
      return
    endif

    nframes = data_size / 2
    allocate(pcm(nframes))
    read(u,iostat=ios) pcm
    close(u)
    if(ios.ne.0) then
      print*,'Failed reading WAV samples: ',trim(fname)
      deallocate(pcm)
      ierr = 5
      return
    endif

    ! Map real 0..4kHz audio to complex baseband at fs=4kHz:
    ! mix by exp(-j*2*pi*2000*t), low-pass filter, then decimate by 3.
    pi = 4.0d0*atan(1.0d0)
    fs = 12000.0d0
    fc = 1800.0d0
    dphi = -2.0d0*pi*(2000.0d0/fs)

    ! 25-tap Hamming-windowed low-pass FIR at 12 kHz.
    hsum = 0.0d0
    do k=-m,m
      if(k.eq.0) then
        hk = 2.0d0*fc/fs
      else
        hk = sin(2.0d0*pi*fc*dble(k)/fs)/(pi*dble(k))
      endif
      win = 0.54d0 - 0.46d0*cos(2.0d0*pi*dble(k+m)/dble(ntap-1))
      h(k+m+1) = hk*win
      hsum = hsum + h(k+m+1)
    enddo
    if(hsum.ne.0.0d0) h = h/hsum

    nout = min(NMAX, nframes/3)
    do i=1,nout
      idx = 1 + (i-1)*3
      if(use_fir) then
        sumr = 0.0d0
        sumi = 0.0d0
        do k=-m,m
          it = idx + k
          if(it.ge.1 .and. it.le.nframes) then
            xr = real(pcm(it))/32768.0
            phi = dphi*dble(it-1)
            xw = h(k+m+1)*dble(xr)
            sumr = sumr + xw*cos(phi)
            sumi = sumi + xw*sin(phi)
          endif
        enddo
        ddout(i) = cmplx(real(sumr),real(sumi))
      else
        xr = real(pcm(idx))/32768.0
        phi = dphi*dble(idx-1)
        ddout(i) = cmplx(xr*cos(real(phi)),xr*sin(real(phi)))
      endif
    enddo
    dialfreq = 2000.0d0

    deallocate(pcm)
  end subroutine read_wav_12k_mono

end program ft8d

