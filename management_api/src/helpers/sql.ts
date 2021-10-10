import { Pool } from "pg";
import { pool } from "../config/pg-pool";
import {
	ColumnForTable,
	InsertableForTable,
	Table,
	UpdatableForTable,
	WhereableForTable,
} from "../zapatos/schema";
import * as db from "../zapatos/src";
import {
	AllType,
	SelectOptionsForTable,
	SelectResultMode,
	sql,
	SQL,
	SQLFragment,
	SQLFragmentsMap,
} from "../zapatos/src";

export async function transactionQuery(
	literals: TemplateStringsArray,
	...expressions: SQL[]
) {
	return await db.transaction(
		pool,
		db.Isolation.ReadCommitted,
		(trasactionPool) => sql(literals, ...expressions).run(trasactionPool),
	);
}

export async function query(
	literals: TemplateStringsArray,
	...expressions: SQL[]
) {
	return await sql(literals, ...expressions).run(pool);
}

export const zap = {
	insert: <T extends Table>(table: T, values: InsertableForTable<T>) =>
		db.insert<T>(table, values).run(pool),

	update: <T extends Table>(
		table: T,
		values: UpdatableForTable<T>,
		where: WhereableForTable<T> | SQLFragment,
	) => db.update<T>(table, values, where).run(pool),

	deletes: <T extends Table>(
		table: T,
		where: WhereableForTable<T> | SQLFragment,
	) => db.deletes<T>(table, where).run(pool),

	select: <
		T extends Table,
		C extends ColumnForTable<T>[] | undefined,
		L extends SQLFragmentsMap | undefined,
		E extends SQLFragmentsMap | undefined,
		M extends SelectResultMode = SelectResultMode.Many
	>(
		table: T,
		where: WhereableForTable<T> | SQLFragment | AllType,
		options?: SelectOptionsForTable<T, C, L, E>,
		mode?: M,
	) => db.select<T, C, L, E, M>(table, where, options, mode).run(pool),

	selectOne: <
		T extends Table,
		C extends ColumnForTable<T>[] | undefined,
		L extends SQLFragmentsMap | undefined,
		E extends SQLFragmentsMap | undefined
	>(
		table: T,
		where: WhereableForTable<T> | SQLFragment | AllType,
		options?: SelectOptionsForTable<T, C, L, E>,
	) => db.selectOne<T, C, L, E>(table, where, options).run(pool),

	sql: (literals: TemplateStringsArray, ...expressions: SQL[]) => {
		return sql(literals, ...expressions).run(pool);
	},

	transaction: {
		container: (
			callback: (
				client: db.TxnClient<db.Isolation.ReadCommitted>,
			) => Promise<any[]>,
		) => {
			return db.transaction(pool, db.Isolation.ReadCommitted, callback);
		},

		sql: (
			transactionPool: Pool,
			literals: TemplateStringsArray,
			...expressions: SQL[]
		) => {
			return sql(literals, ...expressions).run(transactionPool);
		},

		single: (literals: TemplateStringsArray, ...expressions: SQL[]) => {
			return db.transaction(
				pool,
				db.Isolation.ReadCommitted,
				(trasactionPool) => sql(literals, ...expressions).run(trasactionPool),
			);
		},
	},
};
