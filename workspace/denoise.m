function signal_p=denoise(out,signal,nspin,n)
for i=0:1:(nspin-1)
     %ѭ��ƽ��8��
    if i > 0
        s = [signal((n+1-i):n) signal(1:(n-i))];
    else
        s = [signal((-i+1):n) signal(1:(-i))];
    end
    iwh=wden(s,'rigrsure','s','mln',3,'sym4');  %����С��ֵ Ӳ��ֵ
    %��ѭ��ƽ��8��
    if -i > 0
        ds = [iwh((n+1+i):n) iwh(1:(n+i))];
    else
        ds = [iwh((i+1):n) iwh(1:(i))];
    end
    out=out+ds;
end
signal_p=out/nspin;