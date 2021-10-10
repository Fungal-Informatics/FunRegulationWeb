/** Types generated for queries found in "src/sql/elaborated.sql" */
import { PreparedQuery } from '@pgtyped/query';

export type Json = null | boolean | number | string | Json[] | { [key: string]: Json };

/** 'SelectBoardsSharedWithUserSql' parameters type */
export interface ISelectBoardsSharedWithUserSqlParams {
  userId: string | null | void;
}

/** 'SelectBoardsSharedWithUserSql' return type */
export interface ISelectBoardsSharedWithUserSqlResult {
  id: string;
  ownerId: string;
  name: string;
  type: string | null;
  createdAt: Date;
  updatedAt: Date;
  isArchived: boolean;
  sharedWith: Json | null;
}

/** 'SelectBoardsSharedWithUserSql' query type */
export interface ISelectBoardsSharedWithUserSqlQuery {
  params: ISelectBoardsSharedWithUserSqlParams;
  result: ISelectBoardsSharedWithUserSqlResult;
}

const selectBoardsSharedWithUserSqlIR: any = {"name":"selectBoardsSharedWithUserSQL","params":[{"name":"userId","transform":{"type":"scalar"},"codeRefs":{"used":[{"a":406,"b":411,"line":19,"col":19}]}}],"usedParamSet":{"userId":true},"statement":{"body":"select \n    b.*, \n    coalesce(\n        json_agg(\n            json_build_object(\n                'id',u.id,\n                'name',u.name\n            )\n        ) filter (where u.id is not null)\n    ) as \"sharedWith\"\nfrom \n    board_sharings bs \nleft join \n    boards b on b.id = bs.\"boardId\"\nleft join \n    users u on bs.\"userId\" = u.id \nwhere \n    bs.\"userId\" = :userId and \n    (bs.\"isArchived\" is false or bs.\"isArchived\" is null)\ngroup by \n    b.id","loc":{"a":42,"b":493,"line":2,"col":0}}};

/**
 * Query generated from SQL:
 * ```
 * select 
 *     b.*, 
 *     coalesce(
 *         json_agg(
 *             json_build_object(
 *                 'id',u.id,
 *                 'name',u.name
 *             )
 *         ) filter (where u.id is not null)
 *     ) as "sharedWith"
 * from 
 *     board_sharings bs 
 * left join 
 *     boards b on b.id = bs."boardId"
 * left join 
 *     users u on bs."userId" = u.id 
 * where 
 *     bs."userId" = :userId and 
 *     (bs."isArchived" is false or bs."isArchived" is null)
 * group by 
 *     b.id
 * ```
 */
export const selectBoardsSharedWithUserSql = new PreparedQuery<ISelectBoardsSharedWithUserSqlParams,ISelectBoardsSharedWithUserSqlResult>(selectBoardsSharedWithUserSqlIR);


/** 'SelectUserOwnedBoardsSql' parameters type */
export interface ISelectUserOwnedBoardsSqlParams {
  userId: string | null | void;
}

/** 'SelectUserOwnedBoardsSql' return type */
export interface ISelectUserOwnedBoardsSqlResult {
  id: string;
  ownerId: string;
  name: string;
  type: string | null;
  createdAt: Date;
  updatedAt: Date;
  isArchived: boolean;
  sharedWith: Json | null;
}

/** 'SelectUserOwnedBoardsSql' query type */
export interface ISelectUserOwnedBoardsSqlQuery {
  params: ISelectUserOwnedBoardsSqlParams;
  result: ISelectUserOwnedBoardsSqlResult;
}

const selectUserOwnedBoardsSqlIR: any = {"name":"selectUserOwnedBoardsSQL","params":[{"name":"userId","transform":{"type":"scalar"},"codeRefs":{"used":[{"a":899,"b":904,"line":45,"col":17}]}}],"usedParamSet":{"userId":true},"statement":{"body":"select \n    b.*,\n    coalesce(\n        json_agg(\n            json_build_object(\n                'id',u.id,\n                'name',u.name\n            )\n        ) filter (where u.id is not null)\n    ) as \"sharedWith\"\nfrom \n    boards b \nleft join \n    board_sharings bs on bs.\"boardId\" = b.id \nleft join \n    users u on bs.\"userId\" = u.id \nwhere \n    \"ownerId\" = :userId and \n    (bs.\"isArchived\" is false or bs.\"isArchived\" is null)\ngroup by \n    b.id","loc":{"a":537,"b":986,"line":28,"col":0}}};

/**
 * Query generated from SQL:
 * ```
 * select 
 *     b.*,
 *     coalesce(
 *         json_agg(
 *             json_build_object(
 *                 'id',u.id,
 *                 'name',u.name
 *             )
 *         ) filter (where u.id is not null)
 *     ) as "sharedWith"
 * from 
 *     boards b 
 * left join 
 *     board_sharings bs on bs."boardId" = b.id 
 * left join 
 *     users u on bs."userId" = u.id 
 * where 
 *     "ownerId" = :userId and 
 *     (bs."isArchived" is false or bs."isArchived" is null)
 * group by 
 *     b.id
 * ```
 */
export const selectUserOwnedBoardsSql = new PreparedQuery<ISelectUserOwnedBoardsSqlParams,ISelectUserOwnedBoardsSqlResult>(selectUserOwnedBoardsSqlIR);


/** 'SelectBoardSql' parameters type */
export interface ISelectBoardSqlParams {
  boardId: string | null | void;
}

/** 'SelectBoardSql' return type */
export interface ISelectBoardSqlResult {
  id: string;
  ownerId: string;
  name: string;
  type: string | null;
  createdAt: Date;
  updatedAt: Date;
  isArchived: boolean;
  sharedWith: Json | null;
}

/** 'SelectBoardSql' query type */
export interface ISelectBoardSqlQuery {
  params: ISelectBoardSqlParams;
  result: ISelectBoardSqlResult;
}

const selectBoardSqlIR: any = {"name":"selectBoardSQL","params":[{"name":"boardId","transform":{"type":"scalar"},"codeRefs":{"used":[{"a":1377,"b":1383,"line":69,"col":14}]}}],"usedParamSet":{"boardId":true},"statement":{"body":"select \n    b.*,\n    coalesce(\n        json_agg(\n            json_build_object(\n                'id',u.id,\n                'name',u.name\n            )\n        ) filter (where u.id is not null)\n    ) as \"sharedWith\"\nfrom \n    boards b \nleft join \n    board_sharings bs on bs.\"boardId\" = b.id \nleft join \n    users u on bs.\"userId\" = u.id \nwhere \n    b.\"id\" = :boardId and \n    (bs.\"isArchived\" is false or bs.\"isArchived\" is null)\ngroup by \n    b.id","loc":{"a":1018,"b":1465,"line":52,"col":0}}};

/**
 * Query generated from SQL:
 * ```
 * select 
 *     b.*,
 *     coalesce(
 *         json_agg(
 *             json_build_object(
 *                 'id',u.id,
 *                 'name',u.name
 *             )
 *         ) filter (where u.id is not null)
 *     ) as "sharedWith"
 * from 
 *     boards b 
 * left join 
 *     board_sharings bs on bs."boardId" = b.id 
 * left join 
 *     users u on bs."userId" = u.id 
 * where 
 *     b."id" = :boardId and 
 *     (bs."isArchived" is false or bs."isArchived" is null)
 * group by 
 *     b.id
 * ```
 */
export const selectBoardSql = new PreparedQuery<ISelectBoardSqlParams,ISelectBoardSqlResult>(selectBoardSqlIR);


