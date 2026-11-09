	subroutine index(n,array,indx)
c	Indexes array(n) such that array(indx(j)) is in
c	ascending order for j = 1,2,,,n. array is not changed.

	dimension array(n),indx(n)

	do j=1,n
	indx(j) = j
	enddo

	k = n/2+1
	ir = n

10	continue
	if (k.gt.1) then
		k = k-1
		indxt = indx(k)
		q= array(indxt)
	else
		indxt = indx(ir)
		q=array(indxt)
		indx(ir) = indx(1)
		ir=ir-1
	if(ir.eq.1) then
		indx(1)=indxt
		return
	endif
	endif
	i = k
	j=k+k
20	if(j.le.ir) then
		if(j.lt.ir) then
		if(array(indx(j)).lt.array(indx(j+1))) j = j+1
		endif
		if(q.lt.array(indx(j))) then
		indx(i)=indx(j)
		i=j
		j=j+j
		else
		j=ir+1
		endif

	go to 20
	endif
	indx(i)=indxt
	goto 10
	end


		