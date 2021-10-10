import { Pool } from "pg";
import { UmzugStorage } from "umzug";
import * as db from "../zapatos/src";

export class DBMigrationStorage implements UmzugStorage {
	private postgreSqlPool: Pool;
	private tableCreationLock: boolean;
	private numberOfRetrys: number;
	private retryCount: number;
	private intervalReference?: NodeJS.Timeout;

	private async createMigrationsTable() {
		await db.sql`
            CREATE TABLE IF NOT EXISTS "migrations"
            ( "name" TEXT PRIMARY KEY);
        `.run(this.postgreSqlPool);
		this.tableCreationLock = false;
	}

	public waitForTableToBeCreated() {
		if (!this.tableCreationLock) return Promise.resolve();
		console.log("DBMigrationStorage: waiting migraitons table to be ready");
		return new Promise((resolve, reject) => {
			this.intervalReference = setInterval(() => {
				if (this.tableCreationLock) {
					console.log("DBMigrationStorage: migrations table is not ready");
					this.retryCount = this.retryCount + 1;
					if (this.retryCount >= this.numberOfRetrys) reject();
				} else {
					console.log("DBMigrationStorage: migrations table is ready");
					if (this.intervalReference) clearInterval(this.intervalReference);
					resolve();
				}
			}, 500);
		});
	}

	constructor(pool: Pool) {
		this.postgreSqlPool = pool;
		this.tableCreationLock = true;
		this.numberOfRetrys = 20;
		this.retryCount = 0;
		this.createMigrationsTable();
	}

	public async logMigration(name: string) {
		await this.waitForTableToBeCreated();
		db.insert(
			//@ts-ignore
			"migrations",
			{ name },
		).run(this.postgreSqlPool);
	}

	public async unlogMigration(name: string) {
		await this.waitForTableToBeCreated();
		db.deletes(
			//@ts-ignore
			"migrations",
			{ name },
		).run(this.postgreSqlPool);
	}

	public async executed() {
		await this.waitForTableToBeCreated();
		const queryResult = await db
			.select(
				//@ts-ignore
				"migrations",
				{},
			)
			.run(this.postgreSqlPool);
		return queryResult.map((m) => m.name).sort();
	}
}
