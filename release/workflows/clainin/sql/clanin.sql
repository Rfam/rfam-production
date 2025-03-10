select
  cm.clan_acc,
  GROUP_CONCAT(f.rfam_id SEPARATOR "\t")
from clan_membership cm, family f
where
  f.rfam_acc = cm.rfam_acc
group by cm.clan_acc
order by cm.clan_acc
