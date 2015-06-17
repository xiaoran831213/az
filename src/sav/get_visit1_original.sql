select
	img_sn
   -- distinct dsc
	-- distinct wgs_img.sbj
from
    wgs_img
        join
    (select 
        sbj, min(visit) as earliest
    from
        wgs_img
    group by sbj) as scr ON wgs_img.sbj = scr.sbj and wgs_img.visit = scr.earliest
where dsc like 'MP%RAGE%' and type = "Original"
order by wgs_img.sbj
