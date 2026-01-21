/*
select
  cm.clan_acc,
  GROUP_CONCAT(f.rfam_id SEPARATOR "\t")
from clan_membership cm, family f
where
  f.rfam_acc = cm.rfam_acc
group by cm.clan_acc
order by cm.clan_acc
*/

SELECT
  cm.clan_acc,
  GROUP_CONCAT(f.rfam_id ORDER BY f.rfam_id SEPARATOR '\t') AS family_ids
FROM clan_membership cm
INNER JOIN family f ON f.rfam_acc = cm.rfam_acc
GROUP BY cm.clan_acc
ORDER BY cm.clan_acc;


