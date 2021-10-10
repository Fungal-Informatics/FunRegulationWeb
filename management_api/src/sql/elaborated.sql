/* @name selectBoardsSharedWithUserSQL */
select 
    b.*, 
    coalesce(
        json_agg(
            json_build_object(
                'id',u.id,
                'name',u.name
            )
        ) filter (where u.id is not null)
    ) as "sharedWith"
from 
    board_sharings bs 
left join 
    boards b on b.id = bs."boardId"
left join 
    users u on bs."userId" = u.id 
where 
    bs."userId" = :userId and 
    (bs."isArchived" is false or bs."isArchived" is null)
group by 
    b.id;




/* @name selectUserOwnedBoardsSQL */
select 
    b.*,
    coalesce(
        json_agg(
            json_build_object(
                'id',u.id,
                'name',u.name
            )
        ) filter (where u.id is not null)
    ) as "sharedWith"
from 
    boards b 
left join 
    board_sharings bs on bs."boardId" = b.id 
left join 
    users u on bs."userId" = u.id 
where 
    "ownerId" = :userId and 
    (bs."isArchived" is false or bs."isArchived" is null)
group by 
    b.id;


/* @name selectBoardSQL */
select 
    b.*,
    coalesce(
        json_agg(
            json_build_object(
                'id',u.id,
                'name',u.name
            )
        ) filter (where u.id is not null)
    ) as "sharedWith"
from 
    boards b 
left join 
    board_sharings bs on bs."boardId" = b.id 
left join 
    users u on bs."userId" = u.id 
where 
    b."id" = :boardId and 
    (bs."isArchived" is false or bs."isArchived" is null)
group by 
    b.id;