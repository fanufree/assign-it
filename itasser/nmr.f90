! requires pulchra.dat
module nmr
    implicit none
    ! variables for adding backbone atoms
    integer,parameter:: num_stat = 2363
    integer,parameter:: num_stat_pro = 1432 ! stats for proline
    integer, dimension(num_stat,3),save:: nco_stat_bin
    integer, dimension(num_stat_pro,3),save:: nco_stat_pro_bin
    real, dimension(num_stat,8,3),save:: nco_stat
    real, dimension(num_stat_pro,8,3),save:: nco_stat_pro
    real,parameter:: INVALIDANGLE = 999
    real,parameter:: M_PI = 3.14159265358979323846
    real,parameter:: pi180 = 180.0/M_PI
    logical,save:: initStatFlag=.false. ! set to true the first time build_bb is called

    ! variables for NOE restraints
    integer,parameter:: MAXSEQLEN = 1000
    integer,parameter:: NUMCONTYPES = 16

    integer,parameter:: NNINDEX = 1
    integer,parameter:: CACAINDEX = 2
    integer,parameter:: CBCBINDEX = 3
    integer,parameter:: SCSCINDEX = 4

    integer,parameter:: NCAINDEX = 5
    integer,parameter:: NCBINDEX = 6
    integer,parameter:: NSCINDEX = 7
    integer,parameter:: CACBINDEX = 8

    integer,parameter:: CASCINDEX = 9
    integer,parameter:: CBSCINDEX = 10
    integer,parameter:: CANINDEX = 11
    integer,parameter:: CBNINDEX = 12

    integer,parameter:: SCNINDEX = 13
    integer,parameter:: CBCAINDEX = 14
    integer,parameter:: SCCAINDEX = 15
    integer,parameter:: SCCBINDEX = 16

    character(len=4), dimension(NUMCONTYPES), parameter:: TYPENAMES = (/'NN  ','CACA','CBCB','SCSC','NCA ', &
                  'NCB ','NSC ','CACB','CASC','CBSC','CAN ','CBN ','SCN ','CBCA','SCCA','SCCB'/)

    integer,parameter:: NTYPE = 1
    integer,parameter:: CATYPE = 2
    integer,parameter:: CBTYPE = 3
    integer,parameter:: SCTYPE = 4
    integer,parameter,dimension(NUMCONTYPES):: ATOMTYPE1 = (/NTYPE,CATYPE,CBTYPE,SCTYPE, NTYPE,NTYPE,NTYPE,CATYPE, &
                                  CATYPE,CBTYPE,CATYPE,CBTYPE, SCTYPE,CBTYPE,SCTYPE,SCTYPE/) ! first atom type of contact
    integer,parameter,dimension(NUMCONTYPES):: ATOMTYPE2 = (/NTYPE,CATYPE,CBTYPE,SCTYPE, CATYPE,CBTYPE,SCTYPE,CBTYPE, &
                                  SCTYPE,SCTYPE,NTYPE,NTYPE, NTYPE,CATYPE,CATYPE,CBTYPE/) ! second atom type of contact

    integer,parameter:: MAXCONTACTPARTNERS = 100
    integer, dimension(MAXSEQLEN,NUMCONTYPES),save:: numContactPartners = 0 ! [res i][TYPE_INDEX] = num contacts incident to i
    integer, dimension(MAXSEQLEN,NUMCONTYPES,MAXCONTACTPARTNERS),save:: contactPartners = 0 ! [res i][TYPE_INDEX][ith partner] = res j
    integer,parameter:: SCOREINDEX = 1 ! used to index into scoreDistBound
    integer,parameter:: DISTINDEX = 2 ! distances
    real, dimension(MAXSEQLEN,MAXSEQLEN,NUMCONTYPES,2),save:: scoreDistBound = 0! stores the score,distbound; matrix not necessarily symmetric
    integer,parameter:: MAXAMBIG = 200 ! maximum number of ambiguous restraints
    integer,parameter:: MAXAMBIGGROUP = 10 ! maximum number of restraints in an ambiguous assignment group
    integer,parameter:: RES1INDEX = 1 ! used for ambigRestraints
    integer,parameter:: RES2INDEX = 2 ! used for ambigRestraints
    integer,parameter:: TYPEINDEX = 3 ! used for ambigRestraints
    integer, dimension(MAXAMBIG,MAXAMBIGGROUP,3),save:: ambigRestraints = 0 ! res1, res2, typeIndex (e.g. NNINDEX);
    real, dimension(MAXAMBIG,MAXAMBIGGROUP,2),save:: ambigScoreDistBound = 0 ! score,distbound
    integer, dimension(MAXAMBIG),save:: ambigCounts = 0 ! for each ambiguous restraint, the number of contacts in the restraint
    real, dimension(MAXAMBIG),save:: ambigPenalty = 0  ! penalty of ambig restr viol, equal to sum of individual restraints
    integer,save:: numAmbig = 0 ! will be set when read restraint file
    real,parameter:: LINPENALTYBOUND = 1.0 ! 2.0 ! if dist-distRestr > 0 && dist-distRestr < LINPENALTYBOUND, a linear penalty is applied

    ! energy terms to save for weight training (used by energy_tot,ESHORT,EHB)
    real,dimension(13),save:: NMR_ETERMS = 0
    ! number of restraint violations (includes ambig restraints), set by noeRestraintPenalty
    integer,save:: violCount = 0
    ! number of restraints (includes ambig restraints), set by readNOEAssignments and readAmbigNOEAssignments
    integer,save:: numRestr = 0

    ! 1        2        3        4         5        6         7        8      9     10     11     12    13
    ! ENOE,    ESHORT3, ESHORT4, ESHORT4a, ESHORT9, ESHORT10, ESHORT11 EHB1c, EHB1, EHB1a, EHB1b, EHB4, ENOEA
    ! ENOE_wt, er1,     er3,     er4,      er5,     er6,      er7,     eh1c,  eh1,  eh1a,  eh1b,  eh4,  ENOEA_wt


    public:: build_bb, build_bb_ij, build_bb_ij_nocheck, readTalos, getPhiPsi, getPhiPsi_ij, &
             torsionRestraintPenalty, readNOEAssignments, readAmbigNOEAssignments, noeRestraintPenalty, & ! noeRestraintPenaltyAnti, &
             strtok

    contains

! builds the backbone atoms N,C,O from the input CA coordinates
! based on algorithm in pulchra
! OXT is stored in ox(Lch+1), oy(Lch+1), oz(Lch+1)
! seq contains the AA type of each residue (in cas.f). Defined by aa
!        aa = (/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE', &
!c            -1    0     1     2     3     4     5     6
!           'PRO','MET','ASP','ASN','LEU', &
!c            7     8     9    10    11
!           'LYS','GLU','GLN','ARG', &
!c           12    13    14    15
!           'HIS','PHE','TYR','TRP','CYX'/)
!c           16    17    18    19    20
subroutine build_bb(cax,cay,caz,nx,ny,nz,cx,cy,cz,ox,oy,oz,Lch,seq)
    implicit none
    integer,intent(in):: Lch
    real,dimension(Lch),intent(in):: cax, cay, caz
    real,dimension(Lch),intent(out):: nx, ny, nz
    real,dimension(Lch),intent(out):: cx, cy, cz
    real,dimension(Lch+1),intent(out):: ox, oy, oz ! +1 for terminal O
    integer,dimension(Lch+1,3):: RBINS
    integer,dimension(Lch),intent(in):: seq ! pro=7
    real,dimension(8):: cacoordsx, cacoordsy, cacoordsz
    real,dimension(8):: tmpcoordsx, tmpcoordsy, tmpcoordsz
    real,dimension(8):: tmpstatx, tmpstaty, tmpstatz
    real:: x1, y1, z1
    real:: x2, y2, z2
    real:: x3, y3, z3
    real:: x4, y4, z4
    real:: besthit, hit
    integer:: bestpos
    integer:: i, j,bin13_1, bin13_2, bin14
    real:: rmsd, total, maxrms
    real,dimension(Lch+10,3):: X_COORDS

    if (.not. initStatFlag) then
        call init_backbone_bins()
        initStatFlag = .true.
    endif

    do i=1,Lch+10
        do j=1,3
            X_COORDS(i,j) = 0.
        end do
    end do

    call prepare_rbins(cax,cay,caz,Lch,RBINS,X_COORDS)

    total = 0
    maxrms = 0.0
    do i=1,Lch+1
        x1 = X_COORDS(i+3,1)
        y1 = X_COORDS(i+3,2)
        z1 = X_COORDS(i+3,3)

        x2 = X_COORDS(i+4,1)
        y2 = X_COORDS(i+4,2)
        z2 = X_COORDS(i+4,3)

        x3 = X_COORDS(i+5,1)
        y3 = X_COORDS(i+5,2)
        z3 = X_COORDS(i+5,3)

        x4 = X_COORDS(i+6,1)
        y4 = X_COORDS(i+6,2)
        z4 = X_COORDS(i+6,3)

        cacoordsx(1) = x1
        cacoordsy(1) = y1
        cacoordsz(1) = z1

        cacoordsx(2) = x2
        cacoordsy(2) = y2
        cacoordsz(2) = z2

        cacoordsx(3) = x3
        cacoordsy(3) = y3
        cacoordsz(3) = z3

        cacoordsx(4) = x4
        cacoordsy(4) = y4
        cacoordsz(4) = z4

        bin13_1 = RBINS(i,1)
        bin13_2 = RBINS(i,2)
        bin14 = RBINS(i,3)

        if (i .ne. 1 .and. seq(i-1) == 7) then ! prev res is pro
           j=1
           besthit=1000.
           bestpos=1
           do
             hit = abs(nco_stat_pro_bin(j,1)-bin13_1)+abs(nco_stat_pro_bin(j,2)-bin13_2)+0.2*abs(nco_stat_pro_bin(j,3)-bin14)
             if (hit<besthit) then
                besthit=hit
                bestpos=j
             endif
             j = j+1
             if (nco_stat_pro_bin(j,1) < 0 .or. hit <= 1e-3) exit
           end do

           do j=1,4
                tmpstatx(j) = nco_stat_pro(bestpos,j,1)
                tmpstaty(j) = nco_stat_pro(bestpos,j,2)
                tmpstatz(j) = nco_stat_pro(bestpos,j,3)
           end do

           do j=1,8
                tmpcoordsx(j) = nco_stat_pro(bestpos,j,1)
                tmpcoordsy(j) = nco_stat_pro(bestpos,j,2)
                tmpcoordsz(j) = nco_stat_pro(bestpos,j,3)
           end do
        else
            j=1
            besthit=1000.
            bestpos=1
            do
                hit = abs(nco_stat_bin(j,1)-bin13_1)+abs(nco_stat_bin(j,2)-bin13_2)+0.2*abs(nco_stat_bin(j,3)-bin14)
                if (hit<besthit) then
                    besthit=hit
                    bestpos=j
                endif
                j = j+1
                if (nco_stat_bin(j,1) < 0 .or. hit<=1e-3) exit
            end do

            do j=1,4
                tmpstatx(j) = nco_stat(bestpos,j,1)
                tmpstaty(j) = nco_stat(bestpos,j,2)
                tmpstatz(j) = nco_stat(bestpos,j,3)
            end do

            do j=1,8
                tmpcoordsx(j) = nco_stat(bestpos,j,1)
                tmpcoordsy(j) = nco_stat(bestpos,j,2)
                tmpcoordsz(j) = nco_stat(bestpos,j,3)
            end do
        endif

        rmsd=superimpose2(cacoordsx,cacoordsy,cacoordsz,tmpstatx,tmpstaty,tmpstatz,4,tmpcoordsx,tmpcoordsy,tmpcoordsz,8)
        total = total + rmsd
        if (rmsd>maxrms) then
            maxrms=rmsd
        endif

        ! add  C,O,N
        if (i .ne. 1) then
            cx(i-1) = tmpcoordsx(5)
            cy(i-1) = tmpcoordsy(5)
            cz(i-1) = tmpcoordsz(5)
            ox(i-1) = tmpcoordsx(6)
            oy(i-1) = tmpcoordsy(6)
            oz(i-1) = tmpcoordsz(6)
        endif
        if (i .ne. Lch+1) then
            nx(i) = tmpcoordsx(7)
            ny(i) = tmpcoordsy(7)
            nz(i) = tmpcoordsz(7)
        else ! terminal oxygen instead of nitrogen
            ox(i) = tmpcoordsx(7)
            oy(i) = tmpcoordsy(7)
            oz(i) = tmpcoordsz(7)
        endif
    end do ! end for each residue
end subroutine

! build_bb_ij will adjust start and end according to the following rules:
!   Residues outside of start, end can have torsion angle changes when residues in
!      start, end (inclusive) change. Therefore, +/- 2 residues are added to start, end
!   end-start+1 must be at least 5 because this is the minimum number needed for construction
!      If not, start and end are enlarged equally in both directions (unless one end is
!      at the boundary in which case it is enlarged more in one direction
!   The first 2 residues and the last 3 of the sequence is used to build the ends of
!      the chain, so if start is 2, then it is set to 1, and if end is set to Lch-2 or
!      Lch-1, then it is set to Lch
subroutine build_bb_ij(cax,cay,caz,nx,ny,nz,cx,cy,cz,ox,oy,oz,start,end,Lch,seq)
    implicit none
    integer,intent(in):: start,end
    integer,intent(in):: Lch
    real,dimension(Lch),intent(in):: cax, cay, caz
    real,dimension(Lch),intent(inout):: nx, ny, nz
    real,dimension(Lch),intent(inout):: cx, cy, cz
    real,dimension(Lch+1),intent(inout):: ox, oy, oz ! +1 for terminal O
    integer,dimension(Lch),intent(in):: seq ! pro=7
    integer:: s, e, d, temp

    if (.not. initStatFlag) then
        call init_backbone_bins()
        if (end-start+1 == Lch) then
            initStatFlag = .true.
        endif
    endif

    ! enlarge interval
    s = start-2
    e = end+2
    if (s < 1) then
        s = 1
    endif
    if (e > Lch) then
        e = Lch
    endif

    ! make sure at least width 5
    if (e-s < 5) then
        temp = 5-e+s
        temp = ceiling(real(temp)/2.0)
        s = s-temp
        if (s < 1) then
            d = 1-s
            s = 1
            e = e+d
            if (e > Lch) then
                e = Lch
            endif
        endif
        e = e+temp
        if (e > Lch) then
            d = e-Lch
            e = Lch
            s = s-d
            if (s < 1) then
                s = 1
            endif
        endif
    endif

    ! round end points if near chain boundary
    if (s == 2) then
        s = 1
    endif
    if (e == Lch-2 .or. e == Lch-1) then
        e = Lch
    endif

    call build_bb_ij_nocheck(cax,cay,caz,nx,ny,nz,cx,cy,cz,ox,oy,oz,s,e,Lch,seq)
end subroutine

! used by build_bb_ij
! see build_bb_ij for conditions on start, end. The conditions are skipped here
subroutine build_bb_ij_nocheck(cax,cay,caz,nx,ny,nz,cx,cy,cz,ox,oy,oz,start,end,Lch,seq)
    implicit none
    integer,intent(in):: start,end
    integer,intent(in):: Lch
    real,dimension(Lch),intent(in):: cax, cay, caz
    real,dimension(Lch),intent(inout):: nx, ny, nz
    real,dimension(Lch),intent(inout):: cx, cy, cz
    real,dimension(Lch+1),intent(inout):: ox, oy, oz ! +1 for terminal O
    integer,dimension(Lch),intent(in):: seq ! pro=7
    integer,dimension(end-start+2,3):: RBINS
    real,dimension(8):: cacoordsx, cacoordsy, cacoordsz
    real,dimension(8):: tmpcoordsx, tmpcoordsy, tmpcoordsz
    real,dimension(8):: tmpstatx, tmpstaty, tmpstatz
    real:: x1, y1, z1
    real:: x2, y2, z2
    real:: x3, y3, z3
    real:: x4, y4, z4
    real:: besthit, hit
    real:: rmsd
    integer:: bestpos, length, si
    integer:: i, j,bin13_1, bin13_2, bin14
    real,dimension(end-start+11,3):: X_COORDS
    length = end-start+1

    if (.not. initStatFlag) then
        call init_backbone_bins()
        if (end-start+1 == Lch) then
            initStatFlag = .true.
        endif
    endif

    do i=1,length+10
        do j=1,3
            X_COORDS(i,j) = 0.
        end do
    end do

    call prepare_rbins_ij(cax,cay,caz,start,end,Lch,RBINS,X_COORDS)

    do i=1,length+1
        x1 = X_COORDS(i+3,1)
        y1 = X_COORDS(i+3,2)
        z1 = X_COORDS(i+3,3)

        x2 = X_COORDS(i+4,1)
        y2 = X_COORDS(i+4,2)
        z2 = X_COORDS(i+4,3)

        x3 = X_COORDS(i+5,1)
        y3 = X_COORDS(i+5,2)
        z3 = X_COORDS(i+5,3)

        x4 = X_COORDS(i+6,1)
        y4 = X_COORDS(i+6,2)
        z4 = X_COORDS(i+6,3)

        cacoordsx(1) = x1
        cacoordsy(1) = y1
        cacoordsz(1) = z1

        cacoordsx(2) = x2
        cacoordsy(2) = y2
        cacoordsz(2) = z2

        cacoordsx(3) = x3
        cacoordsy(3) = y3
        cacoordsz(3) = z3

        cacoordsx(4) = x4
        cacoordsy(4) = y4
        cacoordsz(4) = z4

        bin13_1 = RBINS(i,1)
        bin13_2 = RBINS(i,2)
        bin14 = RBINS(i,3)

        si = i+start-1 ! index between star, -end

        if (si .ne. 1 .and. seq(si-1) == 7) then ! prev res is pro
           j=1
           besthit=1000.
           bestpos=1
           do
             hit = abs(nco_stat_pro_bin(j,1)-bin13_1)+abs(nco_stat_pro_bin(j,2)-bin13_2)+0.2*abs(nco_stat_pro_bin(j,3)-bin14)
             if (hit<besthit) then
                besthit=hit
                bestpos=j
             endif
             j = j+1
             if (nco_stat_pro_bin(j,1) < 0 .or. hit <= 1e-3) exit
           end do

           do j=1,4
                tmpstatx(j) = nco_stat_pro(bestpos,j,1)
                tmpstaty(j) = nco_stat_pro(bestpos,j,2)
                tmpstatz(j) = nco_stat_pro(bestpos,j,3)
           end do

           do j=1,8
                tmpcoordsx(j) = nco_stat_pro(bestpos,j,1)
                tmpcoordsy(j) = nco_stat_pro(bestpos,j,2)
                tmpcoordsz(j) = nco_stat_pro(bestpos,j,3)
           end do
        else
            j=1
            besthit=1000.
            bestpos=1
            do
                hit = abs(nco_stat_bin(j,1)-bin13_1)+abs(nco_stat_bin(j,2)-bin13_2)+0.2*abs(nco_stat_bin(j,3)-bin14)
                if (hit<besthit) then
                    besthit=hit
                    bestpos=j
                endif
                j = j+1
                if (nco_stat_bin(j,1) < 0 .or. hit<=1e-3) exit
            end do

            do j=1,4
                tmpstatx(j) = nco_stat(bestpos,j,1)
                tmpstaty(j) = nco_stat(bestpos,j,2)
                tmpstatz(j) = nco_stat(bestpos,j,3)
            end do

            do j=1,8
                tmpcoordsx(j) = nco_stat(bestpos,j,1)
                tmpcoordsy(j) = nco_stat(bestpos,j,2)
                tmpcoordsz(j) = nco_stat(bestpos,j,3)
            end do
        endif

        rmsd = superimpose2(cacoordsx,cacoordsy,cacoordsz,tmpstatx,tmpstaty,tmpstatz,4,tmpcoordsx,tmpcoordsy,tmpcoordsz,8)

        ! add  C,O,N
        if (si .ne. 1) then
            cx(si-1) = tmpcoordsx(5)
            cy(si-1) = tmpcoordsy(5)
            cz(si-1) = tmpcoordsz(5)
            ox(si-1) = tmpcoordsx(6)
            oy(si-1) = tmpcoordsy(6)
            oz(si-1) = tmpcoordsz(6)
        endif
        if (si .ne. Lch+1) then
            nx(si) = tmpcoordsx(7)
            ny(si) = tmpcoordsy(7)
            nz(si) = tmpcoordsz(7)
        else ! terminal oxygen instead of nitrogen
            ox(si) = tmpcoordsx(7)
            oy(si) = tmpcoordsy(7)
            oz(si) = tmpcoordsz(7)
        endif
    end do ! end for each residue
end subroutine

! helper function for build_bb
! populates RBINS and X_COORDS based on input ca coordinates
subroutine prepare_rbins(cax,cay,caz,Lch,RBINS,X_COORDS)
    implicit none
    integer,intent(in):: Lch
    real,dimension(Lch),intent(in):: cax, cay, caz
    integer,dimension(Lch+1,3),intent(out):: RBINS
    integer:: i, bin13_1, bin13_2, bin14
    real:: x1, y1, z1
    real:: x2, y2, z2
    real:: x3, y3, z3
    real:: x4, y4, z4
    real:: r13_1, r13_2, r14
    real,dimension(8):: cacoordsx, cacoordsy, cacoordsz
    real,dimension(8):: tmpcoordsx, tmpcoordsy, tmpcoordsz
    real,dimension(8):: tmpstatx, tmpstaty, tmpstatz
    real,dimension(Lch+10,3),intent(inout):: X_COORDS
    real:: rmsd

    do i=1,Lch
        X_COORDS(i+5,1) = cax(i);
        X_COORDS(i+5,2) = cay(i);
        X_COORDS(i+5,3) = caz(i);
    end do

    ! for rebuilding N-term end (build 2 pseudo res before N-term by superposing
    ! res 3-5 to 1-3, and then keeping the shifted res 1-2 as res -1, 0)
    do i=1,5
       tmpcoordsx(i) = X_COORDS(i+5,1)
       tmpcoordsy(i) = X_COORDS(i+5,2)
       tmpcoordsz(i) = X_COORDS(i+5,3)
    end do

    do i=1,3
        cacoordsx(i) = X_COORDS(i+7,1)
        cacoordsy(i) = X_COORDS(i+7,2)
        cacoordsz(i) = X_COORDS(i+7,3)
    end do

    do i=1,3
        tmpstatx(i) = X_COORDS(i+5,1)
        tmpstaty(i) = X_COORDS(i+5,2)
        tmpstatz(i) = X_COORDS(i+5,3)
    end do

    rmsd = superimpose2(tmpstatx,tmpstaty,tmpstatz,cacoordsx,cacoordsy,cacoordsz,3,tmpcoordsx,tmpcoordsy,tmpcoordsz,5)

    do i=1,2
        X_COORDS(i+3,1) = tmpcoordsx(i)
        X_COORDS(i+3,2) = tmpcoordsy(i)
        X_COORDS(i+3,3) = tmpcoordsz(i)
    end do

    ! for rebuilding C-term end (build 3 pseudo res after C-term by superposing
    ! res [L-4,L-2] to [L-2,L] and then keeping the shifted [L-2,L] as res [L+1,L+3])
    do i=1,5
        tmpcoordsx(i) = X_COORDS(i+Lch,1)
        tmpcoordsy(i) = X_COORDS(i+Lch,2)
        tmpcoordsz(i) = X_COORDS(i+Lch,3)
    end do

    do i=1,3
        cacoordsx(i) = X_COORDS(i+Lch,1)
        cacoordsy(i) = X_COORDS(i+Lch,2)
        cacoordsz(i) = X_COORDS(i+Lch,3)
    end do

    do i=1,3
        tmpstatx(i) = X_COORDS(i+Lch+2,1)
        tmpstaty(i) = X_COORDS(i+Lch+2,2)
        tmpstatz(i) = X_COORDS(i+Lch+2,3)
    end do

    rmsd = superimpose2(tmpstatx,tmpstaty,tmpstatz,cacoordsx,cacoordsy,cacoordsz,3,tmpcoordsx,tmpcoordsy,tmpcoordsz,5)

    do i=1,3
        X_COORDS(i+Lch+5,1) = tmpcoordsx(i+3)
        X_COORDS(i+Lch+5,2) = tmpcoordsy(i+3)
        X_COORDS(i+Lch+5,3) = tmpcoordsz(i+3)
    end do

     ! use ca coords of i-2,i-1,i,i+1 to get res i's bin
    do i=1,Lch+1
        x1 = X_COORDS(i+3,1)
        y1 = X_COORDS(i+3,2)
        z1 = X_COORDS(i+3,3)

        x2 = X_COORDS(i+4,1)
        y2 = X_COORDS(i+4,2)
        z2 = X_COORDS(i+4,3)

        x3 = X_COORDS(i+5,1)
        y3 = X_COORDS(i+5,2)
        z3 = X_COORDS(i+5,3)

        x4 = X_COORDS(i+6,1)
        y4 = X_COORDS(i+6,2)
        z4 = X_COORDS(i+6,3)

        r13_1 = calc_distance(x1, y1, z1, x3, y3, z3)
        r13_2 = calc_distance(x2, y2, z2, x4, y4, z4)
        r14 = calc_r14(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

        bin13_1 = int((r13_1-4.6)/0.3)
        bin13_2 = int((r13_2-4.6)/0.3)
        bin14 = int((r14+11.)/0.3)

        if (bin13_1<0) bin13_1=0
        if (bin13_2<0) bin13_2=0
        if (bin14<0) bin14=0
        if (bin13_1>9) bin13_1=9
        if (bin13_2>9) bin13_2=9
        if (bin14>73) bin14=73
        RBINS(i,1) = bin13_1
        RBINS(i,2) = bin13_2
        RBINS(i,3) = bin14
    end do
end subroutine

! helper function for build_bb_ij
subroutine prepare_rbins_ij(cax,cay,caz,start,end,Lch,RBINS,X_COORDS)
    implicit none
    integer,intent(in):: start, end
    integer,intent(in):: Lch
    real,dimension(Lch),intent(in):: cax, cay, caz
    integer,dimension(end-start+2,3),intent(out):: RBINS
    integer:: i, bin13_1, bin13_2, bin14, length
    real:: x1, y1, z1
    real:: x2, y2, z2
    real:: x3, y3, z3
    real:: x4, y4, z4
    real:: r13_1, r13_2, r14
    real,dimension(8):: cacoordsx, cacoordsy, cacoordsz
    real,dimension(8):: tmpcoordsx, tmpcoordsy, tmpcoordsz
    real,dimension(8):: tmpstatx, tmpstaty, tmpstatz
    real,dimension(end-start+11,3),intent(inout):: X_COORDS
    real:: rmsd

    length = end-start+1

    do i=1,length
        X_COORDS(i+5,1) = cax(i+start-1);
        X_COORDS(i+5,2) = cay(i+start-1);
        X_COORDS(i+5,3) = caz(i+start-1);
    end do

    if (start < 3) then
        ! for rebuilding N-term end (build 2 pseudo res before N-term by superposing
        ! res 3-5 to 1-3, and then keeping the shifted res 1-2 as res -1, 0)
        do i=1,5
            tmpcoordsx(i) = X_COORDS(i+5,1)
            tmpcoordsy(i) = X_COORDS(i+5,2)
            tmpcoordsz(i) = X_COORDS(i+5,3)
        end do

        do i=1,3
            cacoordsx(i) = X_COORDS(i+7,1)
            cacoordsy(i) = X_COORDS(i+7,2)
            cacoordsz(i) = X_COORDS(i+7,3)
        end do

        do i=1,3
            tmpstatx(i) = X_COORDS(i+5,1)
            tmpstaty(i) = X_COORDS(i+5,2)
            tmpstatz(i) = X_COORDS(i+5,3)
        end do

        rmsd = superimpose2(tmpstatx,tmpstaty,tmpstatz,cacoordsx,cacoordsy,cacoordsz,3,tmpcoordsx,tmpcoordsy,tmpcoordsz,5)

        do i=1,2
            X_COORDS(i+3,1) = tmpcoordsx(i)
            X_COORDS(i+3,2) = tmpcoordsy(i)
            X_COORDS(i+3,3) = tmpcoordsz(i)
        end do
    else
        do i=1,2
            X_COORDS(i+3,1) = cax(start-3+i)
            X_COORDS(i+3,2) = cay(start-3+i)
            X_COORDS(i+3,3) = caz(start-3+i)
        end do
    endif

    if (end > Lch-3) then
        ! for rebuilding C-term end (build 3 pseudo res after C-term by superposing
        ! res [L-4,L-2] to [L-2,L] and then keeping the shifted [L-2,L] as res [L+1,L+3])
        do i=1,5
            tmpcoordsx(i) = X_COORDS(i+length,1)
            tmpcoordsy(i) = X_COORDS(i+length,2)
            tmpcoordsz(i) = X_COORDS(i+length,3)
        end do

        do i=1,3
            cacoordsx(i) = X_COORDS(i+length,1)
            cacoordsy(i) = X_COORDS(i+length,2)
            cacoordsz(i) = X_COORDS(i+length,3)
        end do

        do i=1,3
            tmpstatx(i) = X_COORDS(i+length+2,1)
            tmpstaty(i) = X_COORDS(i+length+2,2)
            tmpstatz(i) = X_COORDS(i+length+2,3)
        end do

        rmsd = superimpose2(tmpstatx,tmpstaty,tmpstatz,cacoordsx,cacoordsy,cacoordsz,3,tmpcoordsx,tmpcoordsy,tmpcoordsz,5)

        do i=1,3
            X_COORDS(i+length+5,1) = tmpcoordsx(i+3)
            X_COORDS(i+length+5,2) = tmpcoordsy(i+3)
            X_COORDS(i+length+5,3) = tmpcoordsz(i+3)
        end do
    else
        do i=1,3
            X_COORDS(i+length+5,1) = cax(end+i)
            X_COORDS(i+length+5,2) = cay(end+i)
            X_COORDS(i+length+5,3) = caz(end+i)
        end do
    endif

    ! use ca coords of i-2,i-1,i,i+1 to get res i's bin
    do i=1,length+1
        x1 = X_COORDS(i+3,1)
        y1 = X_COORDS(i+3,2)
        z1 = X_COORDS(i+3,3)

        x2 = X_COORDS(i+4,1)
        y2 = X_COORDS(i+4,2)
        z2 = X_COORDS(i+4,3)

        x3 = X_COORDS(i+5,1)
        y3 = X_COORDS(i+5,2)
        z3 = X_COORDS(i+5,3)

        x4 = X_COORDS(i+6,1)
        y4 = X_COORDS(i+6,2)
        z4 = X_COORDS(i+6,3)

        r13_1 = calc_distance(x1, y1, z1, x3, y3, z3)
        r13_2 = calc_distance(x2, y2, z2, x4, y4, z4)
        r14 = calc_r14(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

        bin13_1 = int((r13_1-4.6)/0.3)
        bin13_2 = int((r13_2-4.6)/0.3)
        bin14 = int((r14+11.)/0.3)

        if (bin13_1<0) bin13_1=0
        if (bin13_2<0) bin13_2=0
        if (bin14<0) bin14=0
        if (bin13_1>9) bin13_1=9
        if (bin13_2>9) bin13_2=9
        if (bin14>73) bin14=73
        RBINS(i,1) = bin13_1
        RBINS(i,2) = bin13_2
        RBINS(i,3) = bin14
    end do
end subroutine

! reads noe.tbl
! initializes
! a) numContactPartners
! b) contactPartners
! c) scoreDistBound
! d) numRestr
subroutine readNOEAssignments(filename)
    implicit none
    character(len=*):: filename
    character(len=50):: line
    character(len=4):: type
    integer:: IOstatus, numContacts, i, res1, res2
    real:: score, dist

    open(unit=64,file=filename,status='old')
    read(64,'(A50)',IOSTAT=IOstatus) line
    read(line,'(I6)') numContacts
    ! write(*,*) 'Num NOE Restraints: ',numContacts

    do i=1,numContacts
        read(64,'(A50)',IOSTAT=IOstatus) line
        if (IOstatus.ne.0) then
            write(*,*) 'Not a NOE restraint formatted file'
            return
        endif
        read(line,100) type,res1,res2,score,dist
        ! add both directions
        if (type.eq.'NN') then
            if (numContactPartners(res1,NNINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,NNINDEX) = numContactPartners(res1,NNINDEX)+1
                contactPartners(res1,NNINDEX,numContactPartners(res1,NNINDEX)) = res2
                scoreDistBound(res1,res2,NNINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,NNINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,NNINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,NNINDEX) = numContactPartners(res2,NNINDEX)+1
                contactPartners(res2,NNINDEX,numContactPartners(res2,NNINDEX)) = res1
                scoreDistBound(res2,res1,NNINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,NNINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'CACA') then
            if (numContactPartners(res1,CACAINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,CACAINDEX) = numContactPartners(res1,CACAINDEX)+1
                contactPartners(res1,CACAINDEX,numContactPartners(res1,CACAINDEX)) = res2
                scoreDistBound(res1,res2,CACAINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,CACAINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,CACAINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,CACAINDEX) = numContactPartners(res2,CACAINDEX)+1
                contactPartners(res2,CACAINDEX,numContactPartners(res2,CACAINDEX)) = res1
                scoreDistBound(res2,res1,CACAINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,CACAINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'CBCB') then
            if (numContactPartners(res1,CBCBINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,CBCBINDEX) = numContactPartners(res1,CBCBINDEX)+1
                contactPartners(res1,CBCBINDEX,numContactPartners(res1,CBCBINDEX)) = res2
                scoreDistBound(res1,res2,CBCBINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,CBCBINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,CBCBINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,CBCBINDEX) = numContactPartners(res2,CBCBINDEX)+1
                contactPartners(res2,CBCBINDEX,numContactPartners(res2,CBCBINDEX)) = res1
                scoreDistBound(res2,res1,CBCBINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,CBCBINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'SCSC') then
            if (numContactPartners(res1,SCSCINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,SCSCINDEX) = numContactPartners(res1,SCSCINDEX)+1
                contactPartners(res1,SCSCINDEX,numContactPartners(res1,SCSCINDEX)) = res2
                scoreDistBound(res1,res2,SCSCINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,SCSCINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,SCSCINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,SCSCINDEX) = numContactPartners(res2,SCSCINDEX)+1
                contactPartners(res2,SCSCINDEX,numContactPartners(res2,SCSCINDEX)) = res1
                scoreDistBound(res2,res1,SCSCINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,SCSCINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'NCA') then
            if (numContactPartners(res1,NCAINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,NCAINDEX) = numContactPartners(res1,NCAINDEX)+1
                contactPartners(res1,NCAINDEX,numContactPartners(res1,NCAINDEX)) = res2
                scoreDistBound(res1,res2,NCAINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,NCAINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,CANINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,CANINDEX) = numContactPartners(res2,CANINDEX)+1
                contactPartners(res2,CANINDEX,numContactPartners(res2,CANINDEX)) = res1
                scoreDistBound(res2,res1,CANINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,CANINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'NCB') then
            if (numContactPartners(res1,NCBINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,NCBINDEX) = numContactPartners(res1,NCBINDEX)+1
                contactPartners(res1,NCBINDEX,numContactPartners(res1,NCBINDEX)) = res2
                scoreDistBound(res1,res2,NCBINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,NCBINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,CBNINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,CBNINDEX) = numContactPartners(res2,CBNINDEX)+1
                contactPartners(res2,CBNINDEX,numContactPartners(res2,CBNINDEX)) = res1
                scoreDistBound(res2,res1,CBNINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,CBNINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'NSC') then
            if (numContactPartners(res1,NSCINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,NSCINDEX) = numContactPartners(res1,NSCINDEX)+1
                contactPartners(res1,NSCINDEX,numContactPartners(res1,NSCINDEX)) = res2
                scoreDistBound(res1,res2,NSCINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,NSCINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,SCNINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,SCNINDEX) = numContactPartners(res2,SCNINDEX)+1
                contactPartners(res2,SCNINDEX,numContactPartners(res2,SCNINDEX)) = res1
                scoreDistBound(res2,res1,SCNINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,SCNINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'CACB') then
            if (numContactPartners(res1,CACBINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,CACBINDEX) = numContactPartners(res1,CACBINDEX)+1
                contactPartners(res1,CACBINDEX,numContactPartners(res1,CACBINDEX)) = res2
                scoreDistBound(res1,res2,CACBINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,CACBINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,CBCAINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,CBCAINDEX) = numContactPartners(res2,CBCAINDEX)+1
                contactPartners(res2,CBCAINDEX,numContactPartners(res2,CBCAINDEX)) = res1
                scoreDistBound(res2,res1,CBCAINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,CBCAINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'CASC') then
            if (numContactPartners(res1,CASCINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,CASCINDEX) = numContactPartners(res1,CASCINDEX)+1
                contactPartners(res1,CASCINDEX,numContactPartners(res1,CASCINDEX)) = res2
                scoreDistBound(res1,res2,CASCINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,CASCINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,SCCAINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,SCCAINDEX) = numContactPartners(res2,SCCAINDEX)+1
                contactPartners(res2,SCCAINDEX,numContactPartners(res2,SCCAINDEX)) = res1
                scoreDistBound(res2,res1,SCCAINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,SCCAINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else if (type.eq.'CBSC') then
            if (numContactPartners(res1,CBSCINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res1,CBSCINDEX) = numContactPartners(res1,CBSCINDEX)+1
                contactPartners(res1,CBSCINDEX,numContactPartners(res1,CBSCINDEX)) = res2
                scoreDistBound(res1,res2,CBSCINDEX,SCOREINDEX) = score
                scoreDistBound(res1,res2,CBSCINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
            if (numContactPartners(res2,SCCBINDEX) < MAXCONTACTPARTNERS) then
                numContactPartners(res2,SCCBINDEX) = numContactPartners(res2,SCCBINDEX)+1
                contactPartners(res2,SCCBINDEX,numContactPartners(res2,SCCBINDEX)) = res1
                scoreDistBound(res2,res1,SCCBINDEX,SCOREINDEX) = score
                scoreDistBound(res2,res1,SCCBINDEX,DISTINDEX) = dist
                numRestr = numRestr + 1
            endif
        else
            cycle
        endif
    end do ! end for each contact in file
    close(64)
    100 format(a4,1x,i5,1x,i5,1x,f9.6,1x,f6.2)
end subroutine

! initializes
! a) ambigRestraints
! b) ambigCounts
! c) ambigPenalty
! d) numAmbig
! e) numRestr
subroutine readAmbigNOEAssignments(filename)
    implicit none
    character(len=*):: filename
    character(len=50):: line
    character(len=4):: type
    integer:: IOstatus, numContacts, i, res1, res2, typeI
    real:: score, dist, scoreSum

    logical:: foundFlag
    open(unit=64,file=filename,status='old')

    do
        read(64,'(A50)',IOSTAT=IOstatus) line ! read AMB
        if (IOstatus.ne.0) then
            if (IOStatus.gt.0) then
                write(*,*) 'Error reading ambig NOE file'
                return
            endif
            exit ! end of file
        endif
        if (len(trim(line)).eq.0) then
            cycle
        endif
        read(line,101) type, numContacts
        if ( (numAmbig+1) > MAXAMBIG) then
            exit;
        endif
        numAmbig = numAmbig+1
        ambigCounts(numAmbig) = numContacts
        scoreSum = 0
        do i=1,numContacts
            if (i > MAXAMBIGGROUP) then
                exit
            endif
            read(64,'(A50)',IOSTAT=IOstatus) line
            if (IOstatus.ne.0) then
                write(*,*) 'Error: Invalid ambig.tbl format'
                return
            endif
            read(line,102) type,res1,res2,score,dist
            ambigScoreDistBound(numAmbig,i,SCOREINDEX) = score
            scoreSum = scoreSum+score
            ambigScoreDistBound(numAmbig,i,DISTINDEX) = dist
            ambigRestraints(numAmbig,i,RES1INDEX) = res1
            ambigRestraints(numAmbig,i,RES2INDEX) = res2

            foundFlag = .false.
            do typeI=1,NUMCONTYPES
                if (TYPENAMES(typeI).eq.type) then
                    foundFlag = .true.
                    exit
                endif
            end do
            if (foundFlag) then
                ambigRestraints(numAmbig,i,TYPEINDEX) = typeI
            else
                write(*,*) 'Error: Invalid type in ambig.tbl'
                return
            endif
        end do ! end for each contact in restraint group
        ambigPenalty(numAmbig) = scoreSum
        read(64,'(A50)',IOSTAT=IOstatus) line ! read END
        if (IOstatus.ne.0) then
            write(*,*) 'Error: Invalid ambig.tbl format. Missing END'
            return
        endif
    end do ! end while read file

    close(64)

    numRestr = numRestr + numAmbig

    101 format(a4,1x,i2)
    102 format(a4,1x,i5,1x,i5,1x,f9.6,1x,f6.2)
end subroutine


! reads talos tab
! phi, psi and their errors in dphi, dpsi
! only the 'Good' predictions (according to TALOS) are stored
! otherwise the prediction is INVALIDANGLE
subroutine readTalos(filename, phiPre, psiPre, dphiPre, dpsiPre, Lch)
    implicit none
    character(len=*):: filename
    integer,intent(in):: Lch
    real,dimension(Lch),intent(out):: phiPre, psiPre, dphiPre, dpsiPre
    character(len=70):: line
    character:: resName
    integer:: IOstatus, resNum, num3, cscount,i,offset
    real:: phi, psi, dphi, dpsi, dist, s2
    character(len=4):: classT, temp

    do i=1,Lch
        phiPre(i) = INVALIDANGLE
        psiPre(i) = INVALIDANGLE
        dphiPre(i) = INVALIDANGLE
        dpsiPre(i) = INVALIDANGLE
    end do

    open(unit=71,file=filename,status='old')

    ! skip to FORMAT string, read OFFSET if available
    offset = 0;
    do
        read(71,120,IOSTAT=IOstatus) line
        if (IOstatus.ne.0) then
            write(*,*) 'Not a TALOS formatted file'
            return
        endif
        if (line(1:6).eq.'OFFSET') then
            temp = line(8:)
            temp = trim(temp)
            read(temp,'(i4)') offset
        endif

        if (line(1:6).eq.'FORMAT') then
            exit
        endif
    end do

    do
        read(71,120, IOSTAT=IOstatus) line
        if (IOstatus.ne.0) then
            if (IOStatus.gt.0) then
             write(*,*) 'Error reading file'
             return
            endif
            exit
        endif
        if (len(trim(line)).eq.0) then
            cycle
        endif

        read(line,130) resNum,resName,phi,psi,dphi,dpsi,dist,s2,num3,cscount,classT
        classT = trim(classT)
        if(classT.eq.'Good') then
            phiPre(resNum+offset) = phi
            psiPre(resNum+offset) = psi
            dphiPre(resNum+offset) = dphi
            dpsiPre(resNum+offset) = dpsi
        else
            cycle
        end if
    end do

    close(71)

    120 format(a70)
    130 format(i4,1x,a1,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f5.3,1x,i2,1x,i2,1x,a4)

    return
end subroutine

!real function noeRestraintPenaltyTest(start,end,Lch,cas,sgs)
!    implicit none
!    integer,intent(in)::start
!    integer,intent(in)::end
!    integer,intent(in)::Lch
!    real,dimension(Lch,3),intent(out):: cas
!    real,dimension(Lch,3),intent(out):: sgs
!    real:: energy
!    real:: xa,ya,za,xg,yg,zg,xp,yp,zp
!    real:: dist, cut
!    integer:: i,j,p,numVio
!
!    energy = 0
!    numVio = 0
!    do i=start,end
!        xg = sgs(i,1)
!        yg = sgs(i,2)
!        zg = sgs(i,3)
!        xa = cas(i,1)
!        ya = cas(i,2)
!        za = cas(i,3)
!        do j=1,numContactPartners(i,SCSCINDEX)
!            p = contactPartners(i,SCSCINDEX,j)
!            if (i < p .or. p < start) then ! to prevent duplicate computation
!                xp = sgs(p,1)
!                yp = sgs(p,2)
!                zp = sgs(p,3)
!                dist=sqrt((xg-xp)**2+(yg-yp)**2+(zg-zp)**2)
!
!                cut = scsc(i,p,DISTINDEX)
!
!                if (dist > cut) then
!                    energy = energy + scsc(i,p,SCOREINDEX)
!                    numVio = numVio +1
!                endif
!            endif
!        end do
!        do j=1,numContactPartners(i,SCCAINDEX)
!            p = contactPartners(i,SCCAINDEX,j)
!            if (i < p .or. p < start) then
!                xp = cas(p,1)
!                yp = cas(p,2)
!                zp = cas(p,3)
!                dist=sqrt((xg-xp)**2+(yg-yp)**2+(zg-zp)**2)
!
!                cut = casc(p,i,DISTINDEX)
!
!                if (dist > cut) then
!                    energy = energy + casc(p,i,SCOREINDEX)
!                    numVio = numVio +1
!                endif
!            endif
!        end do
!        do j=1,numContactPartners(i,CASCINDEX)
!            p = contactPartners(i,CASCINDEX,j)
!            if (i < p .or. p < start) then
!                xp = sgs(p,1)
!                yp = sgs(p,2)
!                zp = sgs(p,3)
!                dist=sqrt((xa-xp)**2+(ya-yp)**2+(za-zp)**2)
!
!                cut = casc(i,p,DISTINDEX)
!
!                if (dist > cut) then
!                    energy = energy + casc(i,p,SCOREINDEX)
!                    numVio = numVio +1
!                endif
!            endif
!        end do
!        do j=1,numContactPartners(i,CACAINDEX)
!            p = contactPartners(i,CACAINDEX,j)
!            if (i < p .or. p < start) then
!                xp = cas(p,1)
!                yp = cas(p,2)
!                zp = cas(p,3)
!                dist=sqrt((xa-xp)**2+(ya-yp)**2+(za-zp)**2)
!
!                cut = caca(i,p,DISTINDEX)+1.5
!
!                if (dist > cut) then
!                    energy = energy + caca(i,p,SCOREINDEX)
!                    numVio = numVio +1
!                endif
!            endif
!        end do
!        do j=1,numContactPartners(i,ASCSINDEX)
!            p = contactPartners(i,ASCSINDEX,j)
!            if (i < p .or. p < start) then ! to prevent duplicate computation
!                xp = sgs(p,1)
!                yp = sgs(p,2)
!                zp = sgs(p,3)
!                dist=sqrt((xg-xp)**2+(yg-yp)**2+(zg-zp)**2)
!
!                cut = ascs(i,p,DISTINDEX)
!
!                if (dist < cut) then
!                    energy = energy + ascs(i,p,SCOREINDEX)
!                    numVio = numVio +1
!                endif
!            endif
!        end do
!    end do
!    noeRestraintPenaltyTest = energy
!    return
!end function

real function noeRestraintPenalty(start,endI,printViol)
    implicit none
    integer,intent(in)::start
    integer,intent(in)::endI
    logical,intent(in)::printViol
    real:: energy
    real:: xi,yi,zi,xp,yp,zp ! coords of contact between residue i and p
    real:: axi,ayi,azi,axp,ayp,azp ! ca of i, p; used for CB, SC
    real:: dist, cut, distDiff
    integer:: i,j,p,contactType,numVio
    integer:: iseq,pseq,nv1,nv2,nv1p,nv2p ! used for CB, SC of residue i, p
    integer:: aTypei, aTypep ! atom type of i and p
    integer:: ambigIndex
    logical:: reComputeFlag ! if true, ambiguous restraint is recomputed
    real:: minDistDiff  ! the minimum restraint violation is used to compute ambig restr energy
    integer, parameter:: ndim=1000
    integer, parameter:: nvec=416

    integer seq(ndim),sec(ndim)
    common/seqe/seq,sec
    integer mv(ndim)
    common/chainm/mv

    ! ica for CB,SC; x for all types except N
    integer ica(0:ndim),x(ndim),y(ndim),z(ndim)
    common/chain1/ica,x,y,z

    ! for CA
    real ex(ndim), ey(ndim), ez(ndim)
    common/echain1/ex,ey,ez

    ! for SC
    real egx(ndim),egy(ndim),egz(ndim)
    common/echain2/egx,egy,egz
    real gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
    common/sg/gx,gy,gz

    ! for CB
    real etx(ndim),ety(ndim),etz(ndim)
    common/echain6/etx,ety,etz
    real hx(nvec,nvec,0:19),hy(nvec,nvec,0:19),hz(nvec,nvec,0:19)
    common/cb/hx,hy,hz

    ! for N
    real cax_coord(ndim),cay_coord(ndim),caz_coord(ndim),nx_coord(ndim),ny_coord(ndim),nz_coord(ndim), &
      cx_coord(ndim),cy_coord(ndim),cz_coord(ndim),ox_coord(ndim),oy_coord(ndim),oz_coord(ndim)
    common/nmr_backbone/cax_coord,cay_coord,caz_coord,nx_coord,ny_coord,nz_coord,cx_coord,cy_coord, &
      cz_coord,ox_coord,oy_coord,oz_coord

    energy = 0
    numVio = 0

    do i=start,endI
        iseq=seq(i)
        do contactType=1,NUMCONTYPES
            if (numContactPartners(i,contactType) > 0) then
                ! get the cartesian coords of aTypei, convert from grid if not NTYPE
                aTypei = ATOMTYPE1(contactType)
                aTypep = ATOMTYPE2(contactType)
                if (aTypei.ne.NTYPE) then
                    if (aTypei.eq.CATYPE) then
                        if(mv(i).gt.0)then
                            xi=real(x(i))
                            yi=real(y(i))
                            zi=real(z(i))
                        else
                            xi=ex(i)
                            yi=ey(i)
                            zi=ez(i)
                        endif
                    else if (aTypei == CBTYPE) then
                        if(mv(i).gt.0)then
                            xi=x(i)+HX(ica(i-1),ica(i),iseq)
                            yi=y(i)+HY(ica(i-1),ica(i),iseq)
                            zi=z(i)+HZ(ica(i-1),ica(i),iseq)
                        else
                            xi=etx(i)
                            yi=ety(i)
                            zi=etz(i)
                        endif
                    else if (aTypei == SCTYPE) then
                        if(mv(i).gt.0)then
                            axi=real(x(i))
                            ayi=real(y(i))
                            azi=real(z(i))
                            nv1=ica(i-1)
                            nv2=ica(i)
                            xi=axi+GX(nv1,nv2,iseq)
                            yi=ayi+GY(nv1,nv2,iseq)
                            zi=azi+GZ(nv1,nv2,iseq)
                        else
                            xi=egx(i)
                            yi=egy(i)
                            zi=egz(i)
                        endif
                    else
                        cycle
                    endif
                    ! convert to cartesian
                    xi = 0.87*xi
                    yi = 0.87*yi
                    zi = 0.87*zi
                else
                    ! use CA coordinates instead for efficiency
                    ! xi = nx_coord(i)
                    ! yi = ny_coord(i)
                    ! zi = nz_coord(i)
                    if(mv(i).gt.0)then
                        xi=real(x(i))
                        yi=real(y(i))
                        zi=real(z(i))
                    else
                        xi=ex(i)
                        yi=ey(i)
                        zi=ez(i)
                    endif
                    ! convert to cartesian if using CA (but not for N)
                    xi = 0.87*xi
                    yi = 0.87*yi
                    zi = 0.87*zi
                endif

                do j=1,numContactPartners(i,contactType)
                    p = contactPartners(i,contactType,j)
                    if (i < p .or. p < start) then ! to prevent duplicate computation
                        pseq=seq(p)
                        ! get the cartesian coords of aTypep, convert from grid if not NTYPE
                        if (aTypep.ne.NTYPE) then
                            if (aTypep.eq.CATYPE) then
                                if(mv(p).gt.0)then
                                    xp=real(x(p))
                                    yp=real(y(p))
                                    zp=real(z(p))
                                else
                                    xp=ex(p)
                                    yp=ey(p)
                                    zp=ez(p)
                                endif
                            else if (aTypep == CBTYPE) then
                                if(mv(p).gt.0)then
                                    xp=x(p)+HX(ica(p-1),ica(p),pseq)
                                    yp=y(p)+HY(ica(p-1),ica(p),pseq)
                                    zp=z(p)+HZ(ica(p-1),ica(p),pseq)
                                else
                                    xp=etx(p)
                                    yp=ety(p)
                                    zp=etz(p)
                                endif
                            else if (aTypep == SCTYPE) then
                                if(mv(p).gt.0)then
                                    axp=real(x(p))
                                    ayp=real(y(p))
                                    azp=real(z(p))
                                    nv1p=ica(p-1)
                                    nv2p=ica(p)
                                    xp=axp+GX(nv1p,nv2p,pseq)
                                    yp=ayp+GY(nv1p,nv2p,pseq)
                                    zp=azp+GZ(nv1p,nv2p,pseq)
                                else
                                    xp=egx(p)
                                    yp=egy(p)
                                    zp=egz(p)
                                endif
                            else
                                cycle
                            endif
                            ! convert to cartesian
                            xp = 0.87*xp
                            yp = 0.87*yp
                            zp = 0.87*zp
                        else
                            ! use CA coordinates instead for efficiency
                            ! xp = nx_coord(p)
                            ! yp = ny_coord(p)
                            ! zp = nz_coord(p)
                            if(mv(p).gt.0)then
                                xp=real(x(p))
                                yp=real(y(p))
                                zp=real(z(p))
                            else
                                xp=ex(p)
                                yp=ey(p)
                                zp=ez(p)
                            endif
                            ! convert to cartesian if using CA (but not for N)
                            xp = 0.87*xp
                            yp = 0.87*yp
                            zp = 0.87*zp
                        endif
                        ! compute energy
                        dist=sqrt((xi-xp)**2+(yi-yp)**2+(zi-zp)**2)
                        cut = scoreDistBound(i,p,contactType,DISTINDEX)
                        distDiff = dist-cut
                        if (distDiff > 0) then
                            if (distDiff < LINPENALTYBOUND) then ! linear penalty
                                energy = energy + scoreDistBound(i,p,contactType,SCOREINDEX)*distDiff/LINPENALTYBOUND
                            else
                                energy = energy + scoreDistBound(i,p,contactType,SCOREINDEX)
                                ! write (*,'(i4,i4,a5,3f12.6)') i,p,TYPENAMES(contactType),dist,cut,distDiff
                                numVio = numVio +1
                                if (printViol) then
                                    write (*,*) 'viol: ',i,p
                                endif
                            endif
                        endif
                    endif ! end if i < p < start
                end do ! end for each contact partner j
            endif ! end if num contact partners > 0
        end do ! end for each contact type
    end do ! end for each residue between start end

    ! compute ambiguous restraint penalty
    do ambigIndex=1,numAmbig
        reComputeFlag = .false.
        do j=1,ambigCounts(ambigIndex)
            ! check if restraint contains a residue in start,end
            i = ambigRestraints(ambigIndex,j,RES1INDEX)
            p = ambigRestraints(ambigIndex,j,RES2INDEX)
            if ( (i >= start.and.i <= endI).or.(p >= start.and.p <= endI) ) then
                ! recompute this ambiguous restraint
                reComputeFlag = .true.
                exit
            endif
        end do
        if (reComputeFlag) then
            minDistDiff = 999999.0
            do j=1,ambigCounts(ambigIndex)
                i = ambigRestraints(ambigIndex,j,RES1INDEX)
                p = ambigRestraints(ambigIndex,j,RES2INDEX)
                contactType = ambigRestraints(ambigIndex,j,TYPEINDEX)
                aTypei = ATOMTYPE1(contactType)
                aTypep = ATOMTYPE2(contactType)
                iseq=seq(i)
                pseq=seq(p)
                if (aTypei.ne.NTYPE) then
                    if (aTypei.eq.CATYPE) then
                        if(mv(i).gt.0)then
                            xi=real(x(i))
                            yi=real(y(i))
                            zi=real(z(i))
                        else
                            xi=ex(i)
                            yi=ey(i)
                            zi=ez(i)
                        endif
                    else if (aTypei == CBTYPE) then
                        if(mv(i).gt.0)then
                            xi=x(i)+HX(ica(i-1),ica(i),iseq)
                            yi=y(i)+HY(ica(i-1),ica(i),iseq)
                            zi=z(i)+HZ(ica(i-1),ica(i),iseq)
                        else
                            xi=etx(i)
                            yi=ety(i)
                            zi=etz(i)
                        endif
                    else if (aTypei == SCTYPE) then
                        if(mv(i).gt.0)then
                            axi=real(x(i))
                            ayi=real(y(i))
                            azi=real(z(i))
                            nv1=ica(i-1)
                            nv2=ica(i)
                            xi=axi+GX(nv1,nv2,iseq)
                            yi=ayi+GY(nv1,nv2,iseq)
                            zi=azi+GZ(nv1,nv2,iseq)
                        else
                            xi=egx(i)
                            yi=egy(i)
                            zi=egz(i)
                        endif
                    else
                        cycle
                    endif
                    ! convert to cartesian
                    xi = 0.87*xi
                    yi = 0.87*yi
                    zi = 0.87*zi
                else
                    ! use CA coordinates instead for efficiency
                    ! xi = nx_coord(i)
                    ! yi = ny_coord(i)
                    ! zi = nz_coord(i)
                    if(mv(i).gt.0)then
                        xi=real(x(i))
                        yi=real(y(i))
                        zi=real(z(i))
                    else
                        xi=ex(i)
                        yi=ey(i)
                        zi=ez(i)
                    endif
                    ! convert to cartesian if using CA (but not for N)
                    xi = 0.87*xi
                    yi = 0.87*yi
                    zi = 0.87*zi
                endif ! end if aTypei
                if (aTypep.ne.NTYPE) then
                    if (aTypep.eq.CATYPE) then
                        if(mv(p).gt.0)then
                            xp=real(x(p))
                            yp=real(y(p))
                            zp=real(z(p))
                        else
                            xp=ex(p)
                            yp=ey(p)
                            zp=ez(p)
                        endif
                    else if (aTypep == CBTYPE) then
                        if(mv(p).gt.0)then
                            xp=x(p)+HX(ica(p-1),ica(p),pseq)
                            yp=y(p)+HY(ica(p-1),ica(p),pseq)
                            zp=z(p)+HZ(ica(p-1),ica(p),pseq)
                        else
                            xp=etx(p)
                            yp=ety(p)
                            zp=etz(p)
                        endif
                    else if (aTypep == SCTYPE) then
                        if(mv(p).gt.0)then
                            axp=real(x(p))
                            ayp=real(y(p))
                            azp=real(z(p))
                            nv1p=ica(p-1)
                            nv2p=ica(p)
                            xp=axp+GX(nv1p,nv2p,pseq)
                            yp=ayp+GY(nv1p,nv2p,pseq)
                            zp=azp+GZ(nv1p,nv2p,pseq)
                        else
                            xp=egx(p)
                            yp=egy(p)
                            zp=egz(p)
                        endif
                    else
                        cycle
                    endif
                    ! convert to cartesian
                    xp = 0.87*xp
                    yp = 0.87*yp
                    zp = 0.87*zp
                else
                    ! use CA coordinates instead for efficiency
                    ! xp = nx_coord(p)
                    ! yp = ny_coord(p)
                    ! zp = nz_coord(p)
                    if(mv(p).gt.0)then
                        xp=real(x(p))
                        yp=real(y(p))
                        zp=real(z(p))
                    else
                        xp=ex(p)
                        yp=ey(p)
                        zp=ez(p)
                    endif
                    ! convert to cartesian if using CA (but not for N)
                    xp = 0.87*xp
                    yp = 0.87*yp
                    zp = 0.87*zp
                endif ! end if aTypep
                dist=sqrt((xi-xp)**2+(yi-yp)**2+(zi-zp)**2)
                cut = ambigScoreDistBound(ambigIndex,j,DISTINDEX)
                distDiff = dist-cut
                if (distDiff > 0) then
                    if (distDiff < minDistDiff) then
                        minDistDiff = distDiff
                    endif
                else
                    minDistDiff = 0 ! restraint satisified
                    exit
                endif
            end do ! end for each contact in ambig restraint
            if (minDistDiff > 0) then
                if (minDistDiff < LINPENALTYBOUND) then
                    energy = energy + ambigPenalty(ambigIndex)*minDistDiff/LINPENALTYBOUND
                else
                    energy = energy + ambigPenalty(ambigIndex)
                    numVio = numVio +1
                    if (printViol) then
                        do j=1,ambigCounts(ambigIndex)
                            i = ambigRestraints(ambigIndex,j,RES1INDEX)
                            p = ambigRestraints(ambigIndex,j,RES2INDEX)
                            write (*,*) 'ambig:',ambigIndex,'viol: ',i,p
                        end do
                    endif
                endif
            endif
        endif ! end if reComputeFlag
    end do ! end for each ambig restraint
    noeRestraintPenalty = energy

    violCount = numVio;

  !  write (*,*) 'NumViols: ',numVio
  !  if (endI-start+1 == 102) then
  !      write (*,*) start,endi,energy
  !  endif;
!    if (end-start > 100) then
!        write (*,'(i5,i5,f12.2)') start, end, energy
!    endif
!    write(*,*) 'NOEEnergy',energy
    return
end function

!real function noeRestraintPenaltyAnti(start,endI)
!    implicit none
!    integer,intent(in)::start
!    integer,intent(in)::endI
!    real:: energy
!    real:: axi,ayi,azi,axp,ayp,azp,agxi,agyi,agzi,agxp,agyp,agzp
!    real:: dist, cut
!    integer:: i,j,p,iseq,pseq,nv1,nv2, numVio
!    integer, parameter:: ndim=1000
!    integer, parameter:: nvec=416
!    integer seq(ndim),sec(ndim)
!    common/seqe/seq,sec
!    integer mv(ndim)
!    common/chainm/mv
!    integer ica(0:ndim),x(ndim),y(ndim),z(ndim)
!    common/chain1/ica,x,y,z
!    real ex(ndim), ey(ndim), ez(ndim)
!    common/echain1/ex,ey,ez
!    real egx(ndim),egy(ndim),egz(ndim)
!    common/echain2/egx,egy,egz
!    real gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
!    common/sg/gx,gy,gz
!    energy = 0
!    numVio = 0
!    do i=start,endI
!        iseq=seq(i)
!        if(mv(i).gt.0)then
!            axi=real(x(i))
!            ayi=real(y(i))
!            azi=real(z(i))
!            nv1=ica(i-1)
!            nv2=ica(i)
!            agxi=axi+GX(nv1,nv2,iseq)
!            agyi=ayi+GY(nv1,nv2,iseq)
!            agzi=azi+GZ(nv1,nv2,iseq)
!        else
!            axi=ex(i)
!            ayi=ey(i)
!            azi=ez(i)
!            agxi=egx(i)
!            agyi=egy(i)
!            agzi=egz(i)
!        endif
!        do j=1,numContactPartners(i,ASCSINDEX)
!            p = contactPartners(i,ASCSINDEX,j)
!            if (i < p .or. p < start) then ! to prevent duplicate computation
!                pseq=seq(p)
!                if(mv(p).gt.0)then
!                    agxp=real(x(p))+gx(ica(p-1),ica(p),pseq)
!                    agyp=real(y(p))+gy(ica(p-1),ica(p),pseq)
!                    agzp=real(z(p))+gz(ica(p-1),ica(p),pseq)
!                else
!                    agxp=egx(p)
!                    agyp=egy(p)
!                    agzp=egz(p)
!                endif
!                dist=0.87*sqrt((agxi-agxp)**2+(agyi-agyp)**2+(agzi-agzp)**2) ! in cartesian instead of grid
!                cut = ascs(i,p,DISTINDEX);
!                if (dist < cut) then
!                    energy = energy + ascs(i,p,SCOREINDEX)
!                    numVio = numVio +1
!                endif
!            endif
!        end do
!    end do
!    noeRestraintPenaltyAnti = energy
!  !  if (endI-start+1 == 102) then
!  !      write (*,*) start,endi,energy
!  !  endif;
!!    if (end-start > 100) then
!!        write (*,'(i5,i5,f12.2)') start, end, energy
!!    endif
!!    write(*,*) 'NOEEnergy',energy
!    return
!end function

real function torsionRestraintPenalty(phi,psi,phiPre,psiPre,start,end,Lch)
    implicit none
    integer,intent(in):: Lch
    integer,intent(in):: start
    integer,intent(in):: end
    real,dimension(Lch),intent(in):: phi, psi, phiPre, psiPre
    integer:: i,s,e
    real:: diff,energy
    s = start
    e = end

    ! don't consider the terminal residues
    if (s < 4) then
        s = 4
    endif
    if (e > Lch-4) then
        e = Lch-4
    endif

    energy = 0;

    do i=s,e
        if (phiPre(i) .ne. INVALIDANGLE .and. phi(i) .ne. INVALIDANGLE) then
            diff = abs(phiPre(i)-phi(i)) ! max diff is 180
            if (diff > 180.0) then
                diff = 360.0-diff
            endif
            if (diff > 90.0) then
                energy = energy+2.0*diff/180.0
            endif
        endif
        if (psiPre(i) .ne. INVALIDANGLE .and. psi(i) .ne. INVALIDANGLE) then
            diff = abs(psiPre(i)-psi(i))
            if (diff > 180.0) then
                diff = 360.0-diff
            endif
            if (diff > 90.0) then
                energy = energy+1.5*diff/180.0  ! psi has larger range
            endif
        endif
    end do
    torsionRestraintPenalty = energy
    return
end function

! computes the phi, psi angles and stores them in phi, psi
! residues with no predictions are given INVALIDANGLE as its angles
subroutine getPhiPsi(cax,cay,caz,nx,ny,nz,cx,cy,cz,phi,psi,Lch)
    implicit none
    integer,intent(in):: Lch
    real,dimension(Lch),intent(in):: cax, cay, caz ! ca coords
    real,dimension(Lch),intent(in):: nx, ny, nz ! backbone n coords
    real,dimension(Lch),intent(in):: cx, cy, cz ! backbone c0 coords
    real,dimension(Lch),intent(out):: phi, psi
    real:: n(3), ca(3), c(3), n1(3), c1(3) ! n1 = Ni+1, c1=Ci-1
    integer:: i
    real:: u(3), v(3), w(3), n13(3), n24(3)

    do i=1,Lch
        phi(i) = INVALIDANGLE
        psi(i) = INVALIDANGLE
    end do

    do i=2,Lch-1
        n = (/nx(i),ny(i),nz(i)/)
        ca = (/cax(i),cay(i),caz(i)/)
        c = (/cx(i),cy(i),cz(i)/)

        n1 = (/nx(i+1),ny(i+1),nz(i+1)/)
        c1 =  (/cx(i-1),cy(i-1),cz(i-1)/)

        ! phi
        u = n-c1
        v = ca-n
        w = c-ca
        ! u x v
        n13(1)=u(2)*v(3)-u(3)*v(2)
        n13(2)=u(3)*v(1)-u(1)*v(3)
        n13(3)=u(1)*v(2)-u(2)*v(1)
        n13 = n13/sqrt(sum(n13**2))
        n24(1)=v(2)*w(3)-v(3)*w(2)
        n24(2)=v(3)*w(1)-v(1)*w(3)
        n24(3)=v(1)*w(2)-v(2)*w(1)
        n24 = n24/sqrt(sum(n24**2))
        phi(i) = sign(acos(dot_product(n13,n24)), dot_product(n13,c-ca))
        phi(i) = pi180*phi(i)
        ! psi
        u = ca-n
        v = c-ca
        w = n1-c
        ! u x v
        n13(1)=u(2)*v(3)-u(3)*v(2)
        n13(2)=u(3)*v(1)-u(1)*v(3)
        n13(3)=u(1)*v(2)-u(2)*v(1)
        n13 = n13/sqrt(sum(n13**2))
        n24(1)=v(2)*w(3)-v(3)*w(2)
        n24(2)=v(3)*w(1)-v(1)*w(3)
        n24(3)=v(1)*w(2)-v(2)*w(1)
        n24 = n24/sqrt(sum(n24**2))
        psi(i) = sign(acos(dot_product(n13,n24)), dot_product(n13,n1-c))
        psi(i) = pi180*psi(i)
    end do
end subroutine

! torsion angles from start+1, to end-1 (inclusive)
subroutine getPhiPsi_ij(cax,cay,caz,nx,ny,nz,cx,cy,cz,phi,psi,start,end,Lch)
    implicit none
    integer,intent(in):: Lch
    real,dimension(Lch),intent(in):: cax, cay, caz ! ca coords
    real,dimension(Lch),intent(in):: nx, ny, nz ! backbone n coords
    real,dimension(Lch),intent(in):: cx, cy, cz ! backbone c0 coords
    real,dimension(Lch),intent(out):: phi, psi
    integer,intent(in)::start
    integer,intent(in)::end
    real:: n(3), ca(3), c(3), n1(3), c1(3) ! n1 = Ni+1, c1=Ci-1
    integer:: i
    real:: u(3), v(3), w(3), n13(3), n24(3)

    do i=start,end
        phi(i) = INVALIDANGLE
        psi(i) = INVALIDANGLE
    end do

    do i=start+1,end-1
        n = (/nx(i),ny(i),nz(i)/)
        ca = (/cax(i),cay(i),caz(i)/)
        c = (/cx(i),cy(i),cz(i)/)

        n1 = (/nx(i+1),ny(i+1),nz(i+1)/)
        c1 =  (/cx(i-1),cy(i-1),cz(i-1)/)

        ! phi
        u = n-c1
        v = ca-n
        w = c-ca
        ! u x v
        n13(1)=u(2)*v(3)-u(3)*v(2)
        n13(2)=u(3)*v(1)-u(1)*v(3)
        n13(3)=u(1)*v(2)-u(2)*v(1)
        n13 = n13/sqrt(sum(n13**2))
        n24(1)=v(2)*w(3)-v(3)*w(2)
        n24(2)=v(3)*w(1)-v(1)*w(3)
        n24(3)=v(1)*w(2)-v(2)*w(1)
        n24 = n24/sqrt(sum(n24**2))
        phi(i) = sign(acos(dot_product(n13,n24)), dot_product(n13,c-ca))
        phi(i) = pi180*phi(i)
        ! psi
        u = ca-n
        v = c-ca
        w = n1-c
        ! u x v
        n13(1)=u(2)*v(3)-u(3)*v(2)
        n13(2)=u(3)*v(1)-u(1)*v(3)
        n13(3)=u(1)*v(2)-u(2)*v(1)
        n13 = n13/sqrt(sum(n13**2))
        n24(1)=v(2)*w(3)-v(3)*w(2)
        n24(2)=v(3)*w(1)-v(1)*w(3)
        n24(3)=v(1)*w(2)-v(2)*w(1)
        n24 = n24/sqrt(sum(n24**2))
        psi(i) = sign(acos(dot_product(n13,n24)), dot_product(n13,n1-c))
        psi(i) = pi180*psi(i)
    end do
end subroutine

! superimposition of two sets for coordinates + optional transformation of tpoints
! coords2,tpoints superpose into -> coords1
! (only tpoints modified, coords1,coords2 values stay the same)
! returns rmsd of coords2->coords1
real function superimpose2(coords1x,coords1y,coords1z, &
    coords2x,coords2y,coords2z, &
    npoints, &
    tpointsx, tpointsy, tpointsz, &
    ntpoints)
    implicit none
    integer, intent(in):: npoints
    integer, intent(in):: ntpoints
    real,dimension(npoints), intent(inout):: coords1x
    real,dimension(npoints), intent(inout):: coords1y
    real,dimension(npoints), intent(inout):: coords1z
    real,dimension(npoints), intent(inout):: coords2x
    real,dimension(npoints), intent(inout):: coords2y
    real,dimension(npoints), intent(inout):: coords2z
    real,dimension(ntpoints), intent(inout):: tpointsx
    real,dimension(ntpoints), intent(inout):: tpointsy
    real,dimension(ntpoints), intent(inout):: tpointsz
    real,dimension(3,3):: mat_s, mat_a, mat_b, mat_g, mat_u, tmp_mat
    real:: val, d, alpha, beta, gamma, x, y, z
    real:: cx1, cy1, cz1, cx2, cy2, cz2, tmpx, tmpy, tmpz
    integer:: i,j,k,n

    cx1=0.
    cy1=0.
    cz1=0.
    cx2=0.
    cy2=0.
    cz2=0.

    do i=1,npoints
        cx1 = cx1+coords1x(i)
        cy1 = cy1+coords1y(i)
        cz1 = cz1+coords1z(i)
        cx2 = cx2+coords2x(i)
        cy2 = cy2+coords2y(i)
        cz2 = cz2+coords2z(i)
    end do

     cx1=cx1/real(npoints)
     cy1=cy1/real(npoints)
     cz1=cz1/real(npoints)

     cx2=cx2/real(npoints)
     cy2=cy2/real(npoints)
     cz2=cz2/real(npoints)

    do i=1,npoints
       coords1x(i)=coords1x(i)-cx1
       coords1y(i)=coords1y(i)-cy1
       coords1z(i)=coords1z(i)-cz1
       coords2x(i)=coords2x(i)-cx2
       coords2y(i)=coords2y(i)-cy2
       coords2z(i)=coords2z(i)-cz2
    end do

    do i=1,ntpoints
       tpointsx(i) = tpointsx(i)-cx2
       tpointsy(i) = tpointsy(i)-cy2
       tpointsz(i) = tpointsz(i)-cz2
    end do

    do i=1,3
        do j=1,3
            if (i==j) then
                mat_s(i,j)=1.0
                mat_a(i,j)=1.0
                mat_b(i,j)=1.0
                mat_g(i,j)=1.0
            else
                mat_s(i,j)=0.0
                mat_a(i,j)=0.0
                mat_b(i,j)=0.0
                mat_g(i,j)=0.0
            endif
            mat_u(i,j)=0.
        end do
    end do

    do n=1,npoints

    end do

    do n=1,npoints
      mat_u(1,1)=mat_u(1,1)+coords1x(n)*coords2x(n)
      mat_u(1,2)=mat_u(1,2)+coords1x(n)*coords2y(n)
      mat_u(1,3)=mat_u(1,3)+coords1x(n)*coords2z(n)
      mat_u(2,1)=mat_u(2,1)+coords1y(n)*coords2x(n)
      mat_u(2,2)=mat_u(2,2)+coords1y(n)*coords2y(n)
      mat_u(2,3)=mat_u(2,3)+coords1y(n)*coords2z(n)
      mat_u(3,1)=mat_u(3,1)+coords1z(n)*coords2x(n)
      mat_u(3,2)=mat_u(3,2)+coords1z(n)*coords2y(n)
      mat_u(3,3)=mat_u(3,3)+coords1z(n)*coords2z(n)
    end do

    do i=1,3
        do j=1,3
            tmp_mat(i,j)=0.
        end do
    end do

    do
        d=mat_u(3,2)-mat_u(2,3)
        if (d==0) then
            alpha=0
        else
            alpha=atan(d/(mat_u(2,2)+mat_u(3,3)))
        endif
        if (cos(alpha)*(mat_u(2,2)+mat_u(3,3))+sin(alpha)*(mat_u(3,2)-mat_u(2,3))<0.0) then
            alpha=alpha+M_PI
        endif
        mat_a(2,2)=cos(alpha)
        mat_a(3,3)=cos(alpha)
        mat_a(3,2)=sin(alpha)
        mat_a(2,3)=-mat_a(3,2)
        do i=1,3
            do j=1,3
                do k=1,3
                    tmp_mat(i,j)=tmp_mat(i,j)+mat_u(i,k)*mat_a(j,k)
                end do
            end do
        end do

        do i=1,3
            do j=1,3
              mat_u(i,j)=tmp_mat(i,j)
              tmp_mat(i,j)=0.
            end do
        end do

       do i=1,3
            do j=1,3
                do k=1,3
                    tmp_mat(i,j)=tmp_mat(i,j)+mat_a(i,k)*mat_s(k,j)
                end do
            end do
        end do

        do i=1,3
            do j=1,3
               mat_s(i,j)=tmp_mat(i,j)
               tmp_mat(i,j)=0.
            end do
        end do

        d=mat_u(1,3)-mat_u(3,1)
        if (d==0) then
            beta=0
        else
            beta=atan(d/(mat_u(1,1)+mat_u(3,3)))
        endif

        if (cos(beta)*(mat_u(1,1)+mat_u(3,3))+sin(beta)*(mat_u(1,3)-mat_u(3,1))<0.0) then
            beta=beta+M_PI
        endif

        mat_b(1,1)=cos(beta)
        mat_b(3,3)=cos(beta)
        mat_b(1,3)=sin(beta)
        mat_b(3,1)=-mat_b(1,3)

        do i=1,3
            do j=1,3
                do k=1,3
                    tmp_mat(i,j)=tmp_mat(i,j)+mat_u(i,k)*mat_b(j,k)
                end do
            end do
        end do

        do i=1,3
            do j=1,3
                mat_u(i,j)=tmp_mat(i,j)
                tmp_mat(i,j)=0.
            end do
        end do

        do i=1,3
            do j=1,3
                do k=1,3
                    tmp_mat(i,j)=tmp_mat(i,j)+mat_b(i,k)*mat_s(k,j)
                end do
            end do
        end do

        do i=1,3
            do j=1,3
                mat_s(i,j)=tmp_mat(i,j)
                tmp_mat(i,j)=0.
            end do
        end do
        d=mat_u(2,1)-mat_u(1,2)

        if (d==0) then
            gamma=0
        else
            gamma=atan(d/(mat_u(1,1)+mat_u(2,2)))
        endif

       if (cos(gamma)*(mat_u(1,1)+mat_u(2,2))+sin(gamma)*(mat_u(2,1)-mat_u(1,2))<0.0) then
          gamma=gamma+M_PI
       endif

       mat_g(1,1)=cos(gamma)
       mat_g(2,2)=cos(gamma)
       mat_g(2,1)=sin(gamma)
       mat_g(1,2)=-mat_g(2,1)

        do i=1,3
            do j=1,3
                do k=1,3
                    tmp_mat(i,j)=tmp_mat(i,j)+mat_u(i,k)*mat_g(j,k)
                end do
            end do
        end do

        do i=1,3
            do j=1,3
                mat_u(i,j)=tmp_mat(i,j)
                tmp_mat(i,j)=0.
            end do
        end do

        do i=1,3
            do j=1,3
                do k=1,3
                    tmp_mat(i,j)=tmp_mat(i,j)+mat_g(i,k)*mat_s(k,j)
                end do
            end do
        end do

        do i=1,3
            do j=1,3
                mat_s(i,j)=tmp_mat(i,j)
                tmp_mat(i,j)=0.
            end do
        end do

        val=abs(alpha)+abs(beta)+abs(gamma)
        if (val<=0.001) exit
    end do

    val=0.
   do i=1,npoints
      x=coords2x(i)
      y=coords2y(i)
      z=coords2z(i)
      tmpx=x*mat_s(1,1)+y*mat_s(1,2)+z*mat_s(1,3)
      tmpy=x*mat_s(2,1)+y*mat_s(2,2)+z*mat_s(2,3)
      tmpz=x*mat_s(3,1)+y*mat_s(3,2)+z*mat_s(3,3)
      x=coords1x(i)-tmpx
      y=coords1y(i)-tmpy
      z=coords1z(i)-tmpz
      val=val+x*x+y*y+z*z
   end do

    do i=1,ntpoints
      x=tpointsx(i)
      y=tpointsy(i)
      z=tpointsz(i)
      tpointsx(i)=x*mat_s(1,1)+y*mat_s(1,2)+z*mat_s(1,3)
      tpointsy(i)=x*mat_s(2,1)+y*mat_s(2,2)+z*mat_s(2,3)
      tpointsz(i)=x*mat_s(3,1)+y*mat_s(3,2)+z*mat_s(3,3)
    end do

   do i=1,npoints
      coords1x(i)=coords1x(i)+cx1
      coords1y(i)=coords1y(i)+cy1
      coords1z(i)=coords1z(i)+cz1
      coords2x(i)=coords2x(i)+cx2
      coords2y(i)=coords2y(i)+cy2
      coords2z(i)=coords2z(i)+cz2
   end do

    do i=1,ntpoints
      tpointsx(i)=tpointsx(i)+cx1
      tpointsy(i)=tpointsy(i)+cy1
      tpointsz(i)=tpointsz(i)+cz1
    end do

  superimpose2 = sqrt(val/real(npoints))
  return
end function

! computes the distance
real function calc_distance(x1,y1,z1,x2,y2,z2)
    implicit none
    real,intent(in):: x1,y1,z1,x2,y2,z2
    real:: dx,dy,dz
    real:: dist2
    dx = (x1) - (x2)
    dy = (y1) - (y2)
    dz = (z1) - (z2)
    if (dx .ne. 0 .or. dy .ne. 0 .or. dz .ne. 0 ) then
       dist2 = dx*dx+dy*dy+dz*dz
       calc_distance = sqrt(dist2)
    else
       calc_distance = 0.0
    endif
    return
end function

! r14 chiral distance
real function calc_r14(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
    implicit none
    real,intent(in):: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    real:: r, dx, dy, dz
    real:: vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3
    real:: hand
    dx = x4-x1
    dy = y4-y1
    dz = z4-z1
    r = sqrt(dx*dx+dy*dy+dz*dz)

    vx1=x2-x1
    vy1=y2-y1
    vz1=z2-z1
    vx2=x3-x2
    vy2=y3-y2
    vz2=z3-z2
    vx3=x4-x3
    vy3=y4-y3
    vz3=z4-z3

    hand = (vy1*vz2-vy2*vz1)*vx3+(vz1*vx2-vz2*vx1)*vy3+(vx1*vy2-vx2*vy1)*vz3
    if (hand<0) then
        r=-r
    endif

    calc_r14 = r
    return
end function

! initial call is strtok(string,len(string),startIndex=1,endIndex=1)
! The returned indices startIndex, endIndex are inclusive
! returns true if new token; false otherwise. If false, startIndex and endIndex are undefined
! subsequent calls is strtok(string,len(string),startIndex=previous endIndex+1,endIndex=previous endIndex+1)
! e.g.
!   startIndex = 1
!   endIndex = 1
!   length = len(line)
!   status = strtok(line,length, startIndex, endIndex)
!   if (status) ...
!   startIndex = endIndex+1
!   endIndex = startIndex
!   status = strtok(line,length,startIndex,endIndex)
logical function strtok(strin,strlen,startIndex,endIndex)
    implicit none
    character(len=strlen), intent(in):: strin
    integer, intent(in):: strlen
    integer, intent(inout):: startIndex, endIndex
    character(len=1), parameter:: tab = char(9)
    character(len=1), parameter:: newline = char(10)
    character(len=1), parameter:: ff = char(12)
    character(len=1), parameter:: cr = char(13)
    character(len=1), parameter:: ws(5) = (/' ',tab,newline,cr,ff /)
    integer:: i,j
    logical:: hasDelim, hasChar=.false.

    if (endIndex.gt.strlen .or. startIndex.gt.strlen) then
        strtok = .false.
        return
    endif

    ! advance startIndex until char no longer white space
    do i=startIndex,strlen
        hasDelim = .false.
        do j=1,5
            if (strin(i:i).eq.ws(j)) then
                hasDelim = .true.
                exit
            endif
        end do
        if (.not.hasDelim) then
            exit
        else
            startIndex = startIndex+1
        endif
    end do

    ! might have reached the end
    if (startIndex.gt.strlen) then
        strtok = .false.
        return
    endif

    ! check if char at startIndex is non-delim
    do j=1,5
        if (strin(startIndex:startIndex).ne.ws(j)) then
            hasChar = .true.
            exit
        endif
    end do
    if (.not.haschar) then
        strtok = .false.
        return
    endif

    ! move forward until encounter first white space or reach end of string
    endIndex = startIndex
    do i=startIndex,strlen
        hasDelim = .false.
        do j=1,5
            if (strin(i:i).eq.ws(j)) then
                hasDelim = .true.
                exit
            endif
        end do
        if (hasDelim) then
            exit
        endif
    end do

    if (hasDelim) then
        endIndex = i-1; ! first white space case
    else
        endIndex = i; ! end of string
    endif
    strtok = .true.
    return
end function strtok

! requires pulchra.dat
subroutine init_backbone_bins()
    implicit none
    character(len=25):: line
    integer:: i,j,b1,b2,b3
    real:: r1,r2,r3

    open(unit=87,file='pulchra.dat',status='old')

    do i=1,num_stat
        read(87,'(a25)') line
        read(line,150) b1, b2, b3
        nco_stat_bin(i,1) = b1
        nco_stat_bin(i,2) = b2
        nco_stat_bin(i,3) = b3
        do j=1,8
             read(87,'(a25)') line
             read(line,160) r1,r2,r3
             nco_stat(i,j,1) = r1
             nco_stat(i,j,2) = r2
             nco_stat(i,j,3) = r3
        end do
    end do
    do i=1,num_stat_pro
        read(87,'(a25)') line
        read(line,150) b1, b2, b3
        nco_stat_pro_bin(i,1) = b1
        nco_stat_pro_bin(i,2) = b2
        nco_stat_pro_bin(i,3) = b3
        do j=1,8
            read(87,'(a25)') line
            read(line,160) r1,r2,r3
            nco_stat_pro(i,j,1) = r1
            nco_stat_pro(i,j,2) = r2
            nco_stat_pro(i,j,3) = r3
        end do
    end do

    close(87)

    150 format(3i4)
    160 format(3f8.3)
end subroutine

end module
